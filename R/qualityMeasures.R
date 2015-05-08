#' Quality control measures for binned read counts
#'
#' Calculate various quality control measures on binned read counts.
#'
#' The Shannon entropy is defined as
#' \eqn{S = - sum( n * log(n) )}, where \eqn{n = reads/sum(reads)}.\cr\cr 
#' Spikyness is defined as \eqn{K = sum(abs(diff(reads))) / sum(reads)}.
#' 
#' @param reads A vector of binned read counts.
#' @name qualityControl
#' @author Aaron Taudt
NULL

#' @describeIn qualityControl Calculate the spikyness of a library
qc.spikyness <- function(reads) {
	reads <- as.vector(reads)
	sum.reads <- sum(reads)
	spikyness <- sum(abs(diff(reads))) / sum.reads
	return(spikyness)
}

#' @describeIn qualityControl Calculate the Shannon entropy of a library
qc.entropy <- function(reads) {
	reads <- as.vector(reads)
	total.reads <- sum(reads)
	n <- reads/total.reads
	shannon.entropy <- -sum( n * log(n) , na.rm=TRUE)
	return(shannon.entropy)
}
