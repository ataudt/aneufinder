# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- c("unmappable","monosomy","disomy","trisomy","tetrasomy",paste(5:20,"-somy", sep=""))
state.labels <- c("unmappable","monosomy","disomy","trisomy","tetrasomy","multisomy")
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.aneufinder.univariate <- "aneufinder.univariate.hmm"
state.colors <- c("mappable"="gray68","unmappable"="gray20","monosomy"="gold3","disomy"="springgreen3","trisomy"="orangered1","tetrasomy"="orangered4","x-somy"="purple3")
 
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
