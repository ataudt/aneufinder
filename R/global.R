#' @useDynLib aneufinder
#' @import GenomicRanges
#' @import IRanges
NULL

# =======================================================
# Some global variables that can be used in all functions
# =======================================================
## Class names
class.univariate.hmm <- "aneuHMM"
class.multivariate.hmm <- "aneuMultiHMM"
class.bivariate.hmm <- "aneuBiHMM"
class.hmm.list <- "aneufinder.hmm.list"

## Colors for plotting
state.colors <- c("mapped"="gray68","zero-inflation"="gray90", "nullsomy"="gray90","monosomy"="darkorchid2","disomy"="springgreen2","trisomy"="red3","tetrasomy"="gold2","multisomy"="deepskyblue2","total"="black")
#' aneufinder color scheme
#'
#' Get the color scheme that is used in the \pkg{\link{aneufinder}} plots.
#' @export
stateColors <- function() { return(state.colors) }

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


