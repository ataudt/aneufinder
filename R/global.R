# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- c("nullsomy","null-mixed","monosomy","disomy","trisomy","tetrasomy","multisomy")
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.aneufinder.hmm <- "aneufinder.hmm"
state.colors <- c("mapped"="gray68","nullsomy"="gray20","null-mixed"="gray30","monosomy"="gold3","disomy"="springgreen3","trisomy"="orangered1","tetrasomy"="orangered4","multisomy"="purple3","total"="black")
get.state.labels <- function() { return(state.labels) }
get.state.colors <- function() { return(state.colors[state.labels]) }
 
# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
fsize <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

fprob <- function(mean, variance) {
	return(mean/variance)
}

fmean <- function(size, prob) {
	return(size/prob - size)
}

fvariance <- function(size, prob) {
	return( (size - prob*size) / prob^2 )
}
