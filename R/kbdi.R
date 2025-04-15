kbdi <- function (date, t, p, h, w, wrs = 150, start.date = NULL) { # requires date to compute mean annual precipitation
      d <- as.POSIXlt(date)
      # Data conversion to inches and fahrenheit
      t <- (t * (9 / 5)) + 32
      p <- p * .0393700787402
      h <- h # relative humidity in percentage
      w <- w # wind velocity (km/h)
      wrs <- wrs * .0393700787402
      start <- start.date
      # ----------------------------------------------------------------------
      if (length(t) != length(p)) {
            stop("Input time series of differing lengths")
      }
      m <- matrix(data = c(t, h, p, w), ncol = 4)
      if (any(is.na(m))) {
            warning("Missing values deleted from the input series")
            #             na <- unique(c(which(is.na(t)), which(is.na(p))))
            #             m <- matrix(data = c(t, h, p, w), ncol = 4)
            na <- unique(which(is.na(m), arr.ind = TRUE)[ ,1])
            t <- t[-na] 
            p <- p[-na] 
            d <- d[-na]
            h <- h[-na]
            w <- w[-na]
      }
      # mean annual precip (M) -----------------------------------------------
      yr <- d$year + 1900 
      ind <- which(table(yr) < 360)
      ind1 <- which(is.na(match(yr, as.numeric(names(ind)))))
      M <- mean(tapply(p[ind1], yr[ind1], FUN=sum, na.rm=TRUE))
      # Week index ------------------------------------------------------------
      if(length(d) %% 7 == 0) {
            ap <- c()
      } else {
            ap <- rep(ceiling(length(d) / 7), length(d) %% 7)
      }
      wks <- c(rep(1:(length(d) / 7), each = 7), ap)
      # Detection of the first rainy week-----------------------------------------------
      if (is.null(start)) {
            wp <- tapply(p, wks, sum)
            ind <- which(wp >= wrs)[1]
            if (length(ind) == 0 | is.na(ind)) {
                  stop(paste("'wrs' parameter not reached. Maximum weekly precipitation was", max(tapply(p, wks, sum)) / .0393700787402, "mm"))
            }
            ind2 <- which(wks == unique(wks[ind]))[1]
            start <- which(p[ind2:length(p)] == 0)[1] + ind2 - 1 # first day after the wet spell leading to the wrs parameter
            print(paste("Index initialization on", as.Date(d[start]), "after a wet spell greater than", round(wrs / .0393700787402), "mm in a week"))
      }
      # Fixed start date (snow melt...)
      else {
            start <- which(d == as.POSIXlt(start))
            if (length(start) == 0) {
                  start = 1
            }
            print(paste("Start date of the index set to", as.Date(d[start])))
      }
      # Data filtering, starts after heavy rain----------------------------------------------
      t <- t[start:length(t)]
      p <- p[start:length(p)] 
      d <- d[start:length(d)] 
      h <- h[start:length(h)] 
      w <- w[start:length(w)] 
      ## Wet spell counter --------------------------------------------------------------------------------
      w.spell <- c() # counter of consecutive rainy days
      n <- c()  # counter of days since last rain (N in Noble, 1980)
      net.rain <- rep(NA, length(p)) 
      N <- rep(NA, length(p)) 
      sat <- TRUE # Logical flag to indicate wether the 0.2 precip has been reached (K&B 1968, p.12, Col.3) # Assumption that soil is saturated to init the index
      for (i in 1:length(p)) {
            if(isTRUE(all.equal(0, p[i], tolerance = .001))) {
                  n <- c(i, n)
                  N[i] <- length(n)
                  w.spell <- c()
                  net.rain[i] <- 0
                  sat <- FALSE
            }
            else {
                  n <- c() # reset N counter
                  N[i] <- 0 
                  w.spell <- c(w.spell,i)
                  if (length(w.spell) == 1) {
                        if (p[i] > .2) {
                              net.rain[i] <- p[i] - .2
                              sat <- TRUE
                        }
                        else {
                              net.rain[i] <- 0
                              sat <- FALSE
                        }
                  }
                  if (length(w.spell) > 1) {
                        if (sat == TRUE & sum(p[w.spell]) > .2) {
                              net.rain[i] <- p[i]
                        }
                        if (sat == FALSE & sum(p[w.spell]) > .2) {
                              net.rain[i] <- sum(p[w.spell]) - .2
                              sat <- TRUE
                        }
                        if (sat == FALSE & sum(p[w.spell]) < .2) {
                              net.rain[i] <- 0
                        }
                  } 
            }
      }
      Q <- c(0, rep(NA, (length(t) - 1)))
      dQ <- rep(NA, (length(t))) 
      Ep <- rep(NA, (length(t))) 
      ## KBDI computing-----------------------------------------------------------------------------------------------------------------
      for (i in 2:length(t)) {
            if (net.rain[i] > 0) {
                  Q[i] <- Q[i-1] - (100 * net.rain[i])
            } else {
                  Q[i] <- Q[i-1]
            }
            # Potential evapotranspiration (K&B, 1968)
            # NOTE Alexander (1990) corregidendum on the original equation!!
            Ep[i] <- ((.968 * exp(.0486 * t[i]) - 8.30) * .001) / (1 + (10.88 * exp(-.0441 * M))) 
            dQ[i] <- (800 - Q[i]) * Ep[i]
            Q[i] <- Q[i] + dQ[i]      
      }
      Q[which(Q < 0)] <- 0
      Q <- (Q / 100) / .0393700787402
      ## McArthur's components
      D <- (.191*(Q + 104) * ((N + 1) ^ 1.5)) / ((3.52 * ((N + 1) ^ 1.5)) + (p / .0393700787402) - 1) # McArthur's Drought Factor 
      FDI <- 2 * exp(-.45 + .987*log(D) - .0345 * h + .0338 * ((t - 32)*(5 / 9)) + .0234 * w)  # Forest Fire Danger
      net.rain <- net.rain / .0393700787402 # conversion of net rainfall to mm
      return(cbind.data.frame('date'=as.Date(d, origin = '1970-01-01'), 'net.precip'= net.rain, "KBDI" = Q, 'McArthurDF' = D, "FFDI"=FDI))
}
# End
