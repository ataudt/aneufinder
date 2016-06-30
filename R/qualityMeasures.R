

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

#' @describeIn qualityControl Calculate the spikiness of a library
qc.spikiness <- function(counts) {
	counts <- as.vector(counts)
	sum.counts <- sum(counts)
	spikiness <- sum(abs(diff(counts))) / sum.counts
	return(spikiness)
}

#' @describeIn qualityControl Calculate the Shannon entropy of a library
qc.entropy <- function(counts) {
	counts <- as.vector(counts)
	total.counts <- sum(counts)
	n <- counts/total.counts
	shannon.entropy <- -sum( n * log(n) , na.rm=TRUE)
	return(shannon.entropy)
}

#' @describeIn qualityControl Calculate the Bhattacharyya distance between the '1-somy' and '2-somy' distribution
#' @importFrom stats dnbinom
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
	dist <- -log(sum(sqrt(stats::dnbinom(x, size=distr['1-somy','size'], prob=distr['1-somy','prob']) * stats::dnbinom(x, size=distr['2-somy','size'], prob=distr['2-somy','prob']))))
	return(dist)
}


#' Obtain a data.frame with quality metrics
#' 
#' Obtain a data.frame with quality metrics from a list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' 
#' @param hmms A list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @return A data.frame with columns
#' @author Aaron Taudt
getQC <- function(hmms) {

  ## Helper function
  null2na <- function(x) {
      if (is.null(x)) {
          return(NA)
      } else {
          return(x)
      }
  }
	hmms <- loadHmmsFromFiles(hmms)
	qframe <- list()
	for (i1 in 1:length(hmms)) {
		hmm <- hmms[[i1]]
		qframe[[i1]] <- data.frame(
		                      total.read.count=sum(hmm$bins$counts),
													binsize=width(hmm$bins)[1],
													avg.read.count=mean(hmm$bins$counts),
													spikiness=null2na(hmm$qualityInfo$spikiness),
													entropy=null2na(hmm$qualityInfo$shannon.entropy),
													complexity=null2na(hmm$qualityInfo$complexity[1]),
													loglik=null2na(hmm$convergenceInfo$loglik),
													num.segments=length(hmm$segments),
													bhattacharyya=qc.bhattacharyya(hmm)
													)
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
#' \item Spikiness
#' \item Entropy
#' \item Number of segments
#' \item Bhattacharrya distance
#' \item Loglikelihood
#' }
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param G An integer vector specifying the number of clusters that are compared. See \code{\link[mclust:Mclust]{Mclust}} for details.
#' @param itmax The maximum number of outer and inner iterations for the \code{\link[mclust:Mclust]{Mclust}} function. See \code{\link[mclust:emControl]{emControl}} for details.
#' @param measures The quality measures that are used for the clustering. Supported is any combination of \code{c('spikiness','entropy','num.segments','bhattacharyya','loglik','complexity','avg.read.count','total.read.count','binsize')}. 
#' @param orderBy The quality measure to order the clusters by. Default is \code{'spikiness'}.
#' @param reverseOrder Logical indicating whether the ordering by \code{orderBy} is reversed.
#' @return A \code{list} with the classification, parameters and the \code{\link[mclust]{Mclust}} fit.
#' @author Aaron Taudt
#' @importFrom mclust Mclust emControl mclustBIC
#' @importFrom stats na.omit
#' @export
#'@examples
#'## Get a list of HMMs
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'cl <- clusterByQuality(files)
#'## Plot the clustering and print the parameters
#'plot(cl$Mclust, what='classification')
#'print(cl$parameters)
#'## Select files from the best 2 clusters for further processing
#'best.files <- unlist(cl$classification[1:2])
#'
clusterByQuality <- function(hmms, G=1:9, itmax=c(100,100), measures=c('spikiness','entropy','num.segments','bhattacharyya','loglik','complexity'), orderBy='spikiness', reverseOrder=FALSE) {
	
	hmms <- loadHmmsFromFiles(hmms)
	df <- getQC(hmms)
	df <- df[measures]
	ptm <- startTimedMessage("clustering ...")
	na.mask <- apply(apply(df, 2, is.na), 2, all)
	not.use <- names(which(na.mask))
	if (length(not.use) > 0) {
	  warning("Not using columns ", paste0(not.use, collapse=', '), " because all values are NA.")
	}
	measures <- setdiff(measures, not.use)
	df <- df[!na.mask]
	df <- stats::na.omit(df)
	if (nrow(df)==0) {
	  stop("No values in data.frame.")
	}
	fit <- mclust::Mclust(df, G=G, control=emControl(itmax=itmax))
	stopTimedMessage(ptm)
	params <- fit$parameters$mean
	if (is.null(dim(params))) {
		params <- as.matrix(params)
		colnames(params) <- measures
	} else {
		params <- t(params)
	}
	if (is.null(names(fit$classification))) {
		classification <- list(names(hmms))
	} else {
		classification <- split(names(fit$classification), fit$classification)
	}
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
