# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- factor(c("nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"), levels=c("nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"))
state.distributions <- factor(c('delta','dnbinom','dnbinom','dnbinom','dnbinom','dnbinom'), levels=c('delta','dgeom','dnbinom','dpois','dbinom'))
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.univariate.hmm <- "aneufinder.univariate.hmm"
class.multivariate.hmm <- "aneufinder.multivariate.hmm"
class.hmm.list <- "aneufinder.hmm.list"
state.colors <- c("mapped"="gray68","nullsomy"="gray20","null-mixed"="gray30","monosomy"="gold3","disomy"="springgreen3","trisomy"="orangered1","tetrasomy"="orangered4","multisomy"="purple3","total"="black")
get.state.labels <- function() { return(state.labels) }
get.state.colors <- function() { return(state.colors[state.labels]) }
 
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


