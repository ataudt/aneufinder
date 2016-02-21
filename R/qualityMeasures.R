

#' Quality control measures for binned read counts
#'
#' Calculate various quality control measures on binned read counts.
#'
#' The Shannon entropy is defined as
#' \eqn{S = - sum( n * log(n) )}, where \eqn{n = counts/sum(counts)}.\cr\cr 
#' Spikyness is defined as \eqn{K = sum(abs(diff(counts))) / sum(counts)}.
#' 
#' @param counts A vector of binned read counts.
#' @param hmm An \code{\link{aneuHMM}} object.
#' @return A numeric.
#' @name qualityControl
#' @author Aaron Taudt
NULL

#' @describeIn qualityControl Calculate the spikyness of a library
qc.spikyness <- function(counts) {
	counts <- as.vector(counts)
	sum.counts <- sum(counts)
	spikyness <- sum(abs(diff(counts))) / sum.counts
	return(spikyness)
}

#' @describeIn qualityControl Calculate the Shannon entropy of a library
qc.entropy <- function(counts) {
	counts <- as.vector(counts)
	total.counts <- sum(counts)
	n <- counts/total.counts
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
	x <- 0:max(max(hmm$bins$counts),500)
	dist <- -log(sum(sqrt(dnbinom(x, size=distr['monosomy','size'], prob=distr['monosomy','prob']) * dnbinom(x, size=distr['disomy','size'], prob=distr['disomy','prob']))))
	return(dist)
}


getQC <- function(hmms) {

	hmms <- loadHmmsFromFiles(hmms)
	qframe <- list()
	for (i1 in 1:length(hmms)) {
		hmm <- hmms[[i1]]
		if (!is.null(hmm$segments)) {
			qframe[[i1]] <- data.frame(total.read.count=sum(hmm$bins$counts),
														binsize=width(hmm$bins)[1],
														avg.read.count=mean(hmm$bins$counts),
														spikyness=hmm$qualityInfo$spikyness,
														entropy=hmm$qualityInfo$shannon.entropy,
														complexity=hmm$qualityInfo$complexity,
														loglik=hmm$convergenceInfo$loglik,
														num.segments=length(hmm$segments),
														bhattacharyya=qc.bhattacharyya(hmm))
		}
	}
	names(qframe) <- names(hmms)
	qframe <- do.call(rbind, qframe)
	return(qframe)

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
#' \item Loglikelihood
#' }
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param G An integer vector specifying the number of clusters that are compared. See \code{\link[mclust:Mclust]{Mclust}} for details.
#' @param itmax The maximum number of outer and inner iterations for the \code{\link[mclust:Mclust]{Mclust}} function. See \code{\link[mclust:emControl]{emControl}} for details.
#' @param measures The quality measures that are used for the clustering. Supported is any combination of \code{c('spikyness','entropy','num.segments','bhattacharyya','loglik','complexity','avg.read.count','total.read.count','binsize')}. 
#' @param orderBy The quality measure to order the clusters by. Default is \code{'spikyness'}.
#' @param reverseOrder Logical indicating whether the ordering by \code{orderBy} is reversed.
#' @return A \code{list} with the classification, parameters and the \code{\link[mclust]{Mclust}} fit.
#' @author Aaron Taudt
#' @importFrom mclust Mclust emControl mclustBIC
#' @export
#'@examples
#'## Get a list of HMMs
#'folder <- system.file("extdata/primary-lung/results_univariate", package="aneufinder")
#'files <- list.files(folder, full.names=TRUE)
#'cl <- clusterByQuality(files)
#'## Plot the clustering and print the parameters
#'plot(cl$Mclust, what='classification')
#'print(cl$parameters)
#'## Select files from the best 2 clusters for further processing
#'best.files <- unlist(cl$classification[1:2])
#'
clusterByQuality <- function(hmms, G=1:9, itmax=c(100,100), measures=c('spikyness','entropy','num.segments','bhattacharyya','loglik'), orderBy='spikyness', reverseOrder=FALSE) {
	
	hmms <- loadHmmsFromFiles(hmms)
	df <- getQC(hmms)
	df <- df[measures]
	ptm <- startTimedMessage("clustering ...")
	fit <- mclust::Mclust(df, G=G, control=emControl(itmax=itmax))
	stopTimedMessage(ptm)
	params <- t(fit$parameters$mean)
	classification <- split(names(fit$classification), fit$classification)
	## Reorder clusters
	if (orderBy %in% measures) {
		index <- order(params[,orderBy], decreasing=reverseOrder)
		classification <- classification[index]
		params <- params[index,] # order params last
	}
	names(classification) <- NULL

	cluster <- list(classification=classification, parameters=params, Mclust=fit)
	return(cluster)

}
