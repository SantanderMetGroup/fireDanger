#' @title Fire Weather Index applied to 1D data
#' 
#' @description Implementation of the Canadian Fire Weather Index System to 1D data
#' 
#' @param months Vector of integers corresponding to the months in the data
#' @param Tm Vector of temperature records (deg. Celsius)
#' @param H Vector of relative humidity records (\%)
#' @param r Vector of last 24-h accumulated precipitation (mm)
#' @param W Vector of wind velocity records (Km/h)
#' @param lat Optional. Latitude of the records (in decimal degrees). Default to 46,
#' applying the default parameters of the original FWI System, developed in Canada. See Daylength Adjustment details.
#' @param what Character vector, indicating the components of the FWI system to be returned.
#' Accepted values include any possible subset of the set \{\code{"FFMC"},\code{"DMC"},\code{"DC"}
#' ,\code{"BUI"},\code{"ISI"},\code{"FWI"},\code{"DSR"}\}. Default to \code{"FWI"}.
#' @param init.pars A numeric vector of length 3 with the initialization values for the
#'  FFMC, DMC and DC components, in this order. Default values as proposed by van Wagner (1987).
#' @param spin.up Non-negative integer indicating the number of days considered for FWI spin-up.
#'  Default to 0 (i.e. no spin-up is considered). See the dedicated Section below for further details.
#'  
#' @importFrom stats complete.cases
#'  
#' @return A matrix with the time (days) arranged in rows and the components selected in \code{what}
#' in columns. The attribute \code{colnames} gives the component ordering.
#' 
#' @section Daylength adjustment factors: 
#' By default, the function applies the original FWI daylength adjustment factors for DC and DMC (van Wagner 1987),
#'  although it is possible to adjust them by as a function of latitude through the argument \code{lat}.
#' The reference values used for each latitudinal range are those indicated in p.71 and Tables A3.1 and A3.2 (p69) in
#' Lawson and Armitage (2008).
#' 
#' @section FWI spin-up:
#' FWI is initialized with some values for FFMC, DMC and DC components. This means that the first values of the series are not reliable,
#'  until the index is iterated over several time steps and stabilizes (typically a few days suffice).
#'  Thus, the first index values can be optionally set to \code{NA}. The number of days at the beginning of the series to be set to \code{NA}
#'   is controlled via the \code{spin.up} argument.
#' 
#' @references
#' \itemize{
#' \item Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' \item van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' \item van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' }
#' @author J. Bedia, partially based on the original FORTRAN code by van Wagner and Pickett (1985)
#' 
#' @export
#' 

fwi1D <- function(months, Tm, H, r, W,
                  lat = 46,
                  what = "FWI",
                  init.pars = c(85, 6, 15) ,
                  spin.up = 0) {
      if (any(c(length(Tm), length(H), length(r), length(W)) != length(months))) {
            stop("Input vector lengths differ")
      }
      what <- match.arg(arg = what,
                        choices = c("FFMC", "DMC", "DC", "ISI", "BUI", "FWI", "DSR"),
                        several.ok = TRUE)
      out <- matrix(nrow = length(months), ncol = length(what), dimnames = list(NULL, what))      
      rm.ind <- which(!complete.cases(Tm, H, r, W))
      non.na.ind <- setdiff(1:length(months), rm.ind)
      mes <- months ## (I feel lazy for renaming...)
      if (length(rm.ind) > 0) {
            warning("Missing values were removed from the time series before computation")
            mes <- mes[-rm.ind]
            Tm <- Tm[-rm.ind]
            H <- H[-rm.ind]
            r <- r[-rm.ind]
            W <- W[-rm.ind]
      }
      # Table A3.1 - Lawson and Armitage, p69
      ref.lats <- c(-90, -30, -10, 10, 30, 90) # L&A p71
      loc <- as.character(findInterval(lat, ref.lats))  
      aux <- switch(loc,
                    "1" = list(c(11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10.0, 11.2, 11.8), 
                               c(6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8)),
                    "2" = list(c(10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2),
                               c(6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8)),
                    "3" = list(rep(9, 12), rep(1.4, 12)),
                    "4" = list(c(7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8),
                               c(-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6)),
                    "5" = list(c(6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0),
                               c(-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6))
      )
      Le <- aux[[1]]
      dlf <- aux[[2]]
      aux <- NULL
      if (any(H > 100)) {
            warning("One or more values of humidity above 100% were corrected")
            H[which(H > 100)] <- 100
      }
      if (any(H < 0)) {
            warning("Some negative values of humidity were corrected")
            H[which(H < 0)] <- 0
      }
      if (any(r < 0)) {
            warning("Some negative values of precipitation were corrected")
            r[which(r < 0)] <- 0
      }
      if (any(W < 0)) {
            warning("Some negative values of wind were corrected")
            W[which(W < 0)] <- 0
      }
      f0 <- c(init.pars[1], rep(NA,length(mes))) 
      p0 <- c(init.pars[2], rep(NA,length(mes))) 
      d0 <- c(init.pars[3], rep(NA,length(mes)))
      ISI <- rep(NA, length(mes))
      BUI <- rep(NA, length(mes))
      FWI <- rep(NA, length(mes))
      for (i in 1:length(mes)) {
            m0 <- (147.2 * (101 - f0[i]))/(59.5 + f0[i])
            if (r[i] > 0.5) {
                  rA <- r[i] - 0.5
                  if (m0 <= 150) {
                        mr <- m0 + 42.5 * rA * exp(-100 / (251 - m0)) *  (1 - exp(-6.93 / rA))
                  } else {
                        mr <- m0 + 42.5 * rA * exp(-100 / (251 - m0)) * (1 - exp(-6.93 / rA)) + (0.0015 * (m0 - 150) ^ 2 * (rA ^ (0.5)))
                  }
                  if (mr > 250) {
                        mr <- 250
                  }
                  m0 <- mr
            }
            Ed <- 0.942 * H[i] ^ 0.679 + 11 * exp((H[i] - 100) / 10) + 0.18 * (21.1 - Tm[i]) * (1 - (1 / exp(0.115 * H[i])))
            Ew <- 0.618 * H[i] ^ 0.753 + 10 * exp((H[i] - 100) / 10) + 0.18 * (21.1 - Tm[i]) * (1 - (1 / exp(0.115 * H[i])))
            if (m0 > Ed) {
                  k0 <- 0.424 * (1 - ((H[i] / 100) ^ 1.7)) + 0.0694 * (W[i] ^ 0.5) * (1 - ((H[i] / 100) ^ 8))
                  kd <- k0 * 0.581 * exp(0.0365 * Tm[i])
                  m <- Ed + (m0 - Ed) * (10 ^ (-kd))
            }
            if (m0 < Ed) {
                  if (m0 < Ew) {
                        k1 <- 0.424 * (1 - ((100 - H[i]) / 100) ^ 1.7) + 0.0694 * (W[i] ^ 0.5) * (1 - ((100 - H[i]) / 100) ^ 8)
                        kw <- k1 * (0.581 * (exp(0.0365 * Tm[i])))
                        m <- Ew - ((Ew - m0) * 10 ^ (-kw))
                  }
            }
            if (Ed >= m0 & m0 >= Ew) {
                  m <- m0
            }
            if (m < 0) {
                  m <- 0
            }
            f0[i + 1] <- 59.5 * (250 - m) / (147.2 + m)
            if (Tm[i] < -1.1) {
                  Tm[i] <- -1.1
            }
            K <- 1.894 * (Tm[i] + 1.1) * (100 - H[i]) * Le[mes[i]] * 1e-06
            if (r[i] > 1.5) {
                  re <- (0.92 * r[i]) - 1.27
                  M0 <- 20 + exp(5.6348 - (p0[i] / 43.43))
                  if (p0[i] <= 33) {
                        b <- 100 / (0.5 + (0.3 * p0[i]))
                  } else if (p0[i] > 65) {
                        b <- (6.2 * log(p0[i])) - 17.2
                  } else if (p0[i] > 33 & p0[i] <= 65) {
                        b <- 14 - 1.3 * log(p0[i])
                  }
                  Mr <- M0 + ((1000 * re)/(48.77 + b * re))
                  pr <- 244.72 - 43.43 * log(Mr - 20)
                  if (pr < 0) {
                        pr <- 0
                  }
                  p0[i + 1] <- pr + 100 * K
            } else {
                  p0[i + 1] <- p0[i] + 100 * K
            }
            if (Tm[i] < -2.8) {
                  Tm[i] <- -2.8
            }
            v <- 0.36 * (Tm[i] + 2.8) + dlf[mes[i]]
            if (v < 0) {
                  v <- 0
            }
            if (r[i] > 2.8) {
                  rd <- 0.83 * r[i] - 1.27
                  q0 <- 800 * exp(-d0[i] / 400)
                  qr <- q0 + 3.937 * rd
                  dr <- 400 * log(800 / qr)
                  if (dr < 0) {
                        dr <- 0
                  }
                  d0[i + 1] <- dr + 0.5 * v
            } else {
                  d0[i + 1] <- d0[i] + 0.5 * v
            }
            fW <- exp(0.05039 * W[i])
            fF <- 91.9 * exp(-0.1386 * m) * (1 + ((m ^ 5.31) / (4.93 * 1e+07)))
            ISI[i] <- 0.208 * fW * fF
            if (p0[i + 1] <= 0.4 * d0[i + 1]) {
                  BUI[i] <- (0.8 * p0[i + 1] * d0[i + 1]) / (p0[i + 1] + 0.4 * d0[i + 1])
            } else if (p0[i + 1] > 0.4 * d0[i + 1]) {
                  BUI[i] <- p0[i + 1] - (1 - (0.8 * d0[i + 1]) / (p0[i + 1] + 0.4 * d0[i + 1])) * (0.92 + (0.0114 * p0[i + 1]) ^ 1.7)
            }
            if (!is.finite(BUI[i]) | BUI[i] < 0) {
                  BUI[i] <- 0
            } else if (BUI[i] > 80) {
                  fD <- 1000/(25 + 108.64 * exp(-0.023 * BUI[i]))
            } else {
                  fD <- (0.626 * BUI[i] ^ 0.809) + 2
            }
            B <- 0.1 * ISI[i] * fD
            if (B > 1) {
                  Slog <- 2.72 * (0.434 * log(B)) ^ 0.647
                  FWI[i] <- exp(Slog)
            } else {
                  FWI[i] <- B
            }
      }
      FFMC <- f0[-1]
      DMC <- p0[-1]
      DC <- d0[-1]
      if ("DSR" %in% what) DSR <- 0.0272 * FWI ^ 1.77
      aux <- vapply(X = colnames(out),
                    FUN = function(x) cbind(get(x)),
                    FUN.VALUE = numeric(length(mes)))
      if (length(non.na.ind) > 0) {
            out[non.na.ind,] <- aux
      } else {
            out <- aux
      }
      aux <- NULL
      if (spin.up > 0) out[1:spin.up,] <- NA
      return(out)
}
# End
