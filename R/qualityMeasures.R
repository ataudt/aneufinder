#' Quality control measures for binned read counts
#'
#' Calculate various quality control measures on binned read counts.
#'
#' The Shannon entropy is defined as
#' \eqn{S = - sum( n * log(n) )}, where \eqn{n = reads/sum(reads)}.\cr\cr 
#' Spikyness is defined as \eqn{K = sum(abs(diff(reads))) / sum(reads)}.
#' 
#' @param reads A vector of binned read counts.
#' @param hmm An \code{\link{aneuHMM}} object.
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

#' @describeIn qualityControl Calculate the Bhattacharyya distance between the 'monosomy' and 'disomy' distribution
qc.bhattacharyya <- function(hmm) {
	if (class(hmm)=='aneuHMM') {
		distr <- hmm$distributions
	} else if (class(hmm)=='aneuBiHMM') {
		distr <- hmm$distributions$minus
	}
  if (is.null(distr)) {
    return(0)
  }
	x <- 0:max(max(hmm$bins$reads),500)
	dist <- -log(sum(sqrt(dnbinom(x, size=distr['monosomy','size'], prob=distr['monosomy','prob']) * dnbinom(x, size=distr['disomy','size'], prob=distr['disomy','prob']))))
	return(dist)
}


quality <- function(hmm) {

	if (!is.null(hmm$segments)) {
    qframe <- data.frame(total.read.count=sum(hmm$bins$reads),
  												binsize=width(hmm$bins)[1],
  												avg.read.count=mean(hmm$bins$reads),
  												spikyness=hmm$qualityInfo$spikyness,
  												entropy=hmm$qualityInfo$shannon.entropy,
  												complexity=hmm$qualityInfo$complexity,
  												loglik=hmm$convergenceInfo$loglik,
  												num.segments=length(hmm$segments),
  												bhattacharyya=qc.bhattacharyya(hmm))
  	return(qframe)
	} else {
    return(NULL)
	}

}
												
#' @importFrom mclust Mclust
clusterByQuality <- function(hmms, G=NULL) {
	
	hmms <- loadHmmsFromFiles(hmms)
	df <- do.call(rbind, lapply(hmms, quality))
	df <- df[c('spikyness','entropy','loglik','num.segments','bhattacharyya')]
  if (is.null(G)) {
	  fit <- Mclust(df)
  } else {
    fit <- Mclust(df, G=G)
  }
	return(fit)

}
