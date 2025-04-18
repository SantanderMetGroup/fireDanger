angstronIndex <-
function (t, rh) {
	if (length(t) != length(r)) {
		stop("Input time series of differing lengths")
	}
      AI <- rep(NA, length(t))
	if (any(is.na(t) | is.na(r))) {
		na <- unique(c(which(is.na(t)), which(is.na(r))))
		t <- t[-na]
		r <- r[-na]
	}
	vec <- rep(NA, length(t))
	for (i in 1:length(t)) {
		vec[i] <- (r[i] / 20) + ((27 - t[i])/ 10)
	}
      AI[-na] <- vec
	return(AI)
}
