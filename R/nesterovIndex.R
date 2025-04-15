nesterovIndex <-
function (t, rh, p, modified = FALSE) {
	t <- t
	r <- rh
	p <- p
	modified <- modified
	p[which(p < 0.1)] <- 0
	if (length(t) != length(r) | length(t) != length(p)) {
		stop("Input time series of differing lengths")
	}
      vec <- rep(NA, length(t))
	if (any(is.na(t) | is.na(r) | is.na(p))) {
		unique(c(which(is.na(t)), which(is.na(r)), which(is.na(r)))) -> na
		t <- t[-na] 
		r <- r[-na] 
		p <- p[-na] 
	}
	rep(NA,length(t)) -> td
	for (i in 1:length(t)) {
		td[i] <- 237.7 * ((17.271 * t[i] / (237.7 + t[i])) + log(r[i]/100)) / (17.271 - ((17.271 * t[i] / (237.7 + t[i])) + log(r[i] / 100)))
	}
	NI <- c(0, rep(NA,(length(t)-1)))
	for (i in 2:length(t)) {
		if (p[i] < 3) {
			NI[i] <- t[i]*(t[i] - td[i]) + NI[i-1]
		}
		if (p[i] >= 3 | t[i] < 0) {
			NI[i] <- 0 
		}
	}
	if (isTRUE(modified)) {
		k <- rep(NA,length(p))
		z <- c(1,.8,.6,.4,.2,.1,0)
		p.int <- cut(p, breaks = c(-Inf, 0, .9, 2.9, 5.9, 14.9, 19, Inf), ordered_result = TRUE)
		for (i in 1:length(p)) {
			k[i] <- z[which(levels(p.int) == as.character(p.int[i]))]
		}
		NI <- NI * k
	}
      if (exists("na")) {
            vec[-na] <- NI
            return(vec)
      } else {
            return(NI)
      }
}
