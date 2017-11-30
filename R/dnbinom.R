
# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
dnbinom.size <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

dnbinom.prob <- function(mean, variance) {
	return(mean/variance)
}

dnbinom.mean <- function(size, prob) {
	return(size/prob - size)
}

dnbinom.variance <- function(size, prob) {
	return( (size - prob*size) / prob^2 )
}

dgeom.prob <- function(mean) {
  return( 1/(1+mean) )
}

dgeom.mean <- function(prob) {
	return( (1-prob)/prob )
}

dgeom.variance <- function(prob) {
	return( (1-prob)/prob^2 )
}

dbinom.size <- function(mean, variance) {
	return( mean^2/(mean-variance) )
}

dbinom.prob <- function(mean, variance) {
	return( (mean-variance)/mean )
}

dbinom.mean <- function(size, prob) {
	return( size*prob )
}

dbinom.variance <- function(size, prob) {
	return( size*prob * (1-prob) )
}


