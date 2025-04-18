fmi <-
function (t, rh) {
	t <- t
	r <- rh
	if (length(t) != length(r)) {
		stop("Input time series of differing lengths")
	}
	fmi <- 10 - 0.25 * (t - r)
	return(fmi)
}
