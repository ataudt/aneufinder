# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- c("unmappable","monosomy","disomy","trisomy","tetrasomy",paste(5:20,"-somy", sep=""))
state.labels <- c("unmappable","monosomy","disomy","trisomy","tetrasomy","multisomy")
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.aneufinder.univariate <- "aneufinder.univariate.hmm"
gcolors <- c("unmappable"="gray68","monosomy"="red1","disomy"="springgreen3","trisomy"="red2","tetrasomy"="red3","x-somy"="red4")
 
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
