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


getQC <- function(hmm) {

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
												
#' Cluster based on quality variables
#'
#' This function uses the \pkg{\link{mclust}} package to cluster the input samples based on various quality measures.
#'
#' The employed quality measures are:
#' \itemize{
#' \item Spikyness
#' \item Entropy
#' \item Number of segments
#' \item Bhattacharrya distance
#' }
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param G An integer vector specifying the number of clusters that are compared. See \code{\link[mclust:Mclust]{Mclust}} for details.
#' @param itmax The maximum number of outer and inner iterations for the \code{\link[mclust:Mclust]{Mclust}} function. See \code{\link[mclust:emControl]{emControl}} for details.
#' @param measures The quality measures that are used for the clustering. Supported is any combination of \code{c('spikyness','entropy','num.segments','bhattacharyya','loglik','complexity','avg.read.count','total.read.count','binsize')}. 
#' @param orderBy The quality measure to order the clusters by. Default is \code{'spikyness'}.
#' @param reverseOrder Logical indicating whether the ordering by \code{orderBy} is reversed.
#' @author Aaron Taudt
#' @importFrom mclust Mclust emControl
#' @export
clusterByQuality <- function(hmms, G=1:9, itmax=c(100,100), measures=c('spikyness','entropy','num.segments','bhattacharyya'), orderBy='spikyness', reverseOrder=FALSE) {
	
	hmms <- loadHmmsFromFiles(hmms)
	df <- do.call(rbind, lapply(hmms, getQC))
	df <- df[measures]
	message("clustering ...", appendLF=F); ptm <- proc.time()
	fit <- Mclust(df, G=G, control=emControl(itmax=itmax))
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	params <- t(fit$parameters$mean)
	classification <- split(names(fit$classification), fit$classification)
	## Reorder clusters
	if (orderBy %in% measures) {
		index <- order(params[,orderBy], decreasing=reverseOrder)
		classification <- classification[index]
		params <- params[index,] # order params last
	}
	names(classification) <- NULL

	cluster <- list(classification=classification, parameters=params, fit=fit)
	return(cluster)

}
