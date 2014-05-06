fsize <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

fprob <- function(mean, variance) {
	return(mean/variance)
}
