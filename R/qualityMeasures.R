

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
  if (is.null(counts)) {
      return(NA)
  }
	counts <- as.vector(counts)
	sum.counts <- sum(counts)
	spikiness <- sum(abs(diff(counts))) / sum.counts
	return(spikiness)
}

#' @describeIn qualityControl Calculate the Shannon entropy of a library
qc.entropy <- function(counts) {
  if (is.null(counts)) {
      return(NA)
  }
	counts <- as.vector(counts)
	total.counts <- sum(counts)
	n <- counts/total.counts
	entropy <- -sum( n * log(n) , na.rm=TRUE)
	return(entropy)
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
    return(NA)
  }
	x <- 0:max(max(hmm$bins$counts),500)
	dist <- -log(sum(sqrt(stats::dnbinom(x, size=distr['1-somy','size'], prob=distr['1-somy','prob']) * stats::dnbinom(x, size=distr['2-somy','size'], prob=distr['2-somy','prob']))))
	return(dist)
}

#' @describeIn qualityControl Sum-of-squares distance from the read counts to the fitted distributions
#' @importFrom stats dnbinom
qc.sos <- function(hmm) {
	if (class(hmm)=='aneuHMM') {
		distr <- hmm$distributions
    mu <- distr$mu
    names(mu) <- rownames(distr)
    sos <- sum( (hmm$bins$counts - mu[as.character(hmm$bins$state)]) ^ 2 )
	} else if (class(hmm)=='aneuBiHMM') {
		distr <- hmm$distributions$minus
    mu <- distr$mu
    names(mu) <- rownames(distr)
    sos <- sum( (c(hmm$bins$mcounts - mu[as.character(hmm$bins$mstate)],
                   hmm$bins$pcounts - mu[as.character(hmm$bins$pstate)])
                 ) ^ 2 )
	}
  if (is.null(distr)) {
    return(NA)
  }
	return(sos)
}


#' Obtain a data.frame with quality metrics
#' 
#' Obtain a data.frame with quality metrics from a list of \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' 
#' The employed quality measures are:
#' \itemize{
#' \item total.read.count: Total read count.
#' \item avg.binsize: Average binsize.
#' \item avg.read.count: Average read count.
#' \item spikiness: Bin-to-bin variability of read count.
#' \item entropy: Shannon entropy of read counts.
#' \item complexity: Library complexity approximated with a Michaelis-Menten curve.
#' \item loglik: Loglikelihood of the Hidden Markov Model.
#' \item num.segments: Number of copy number segments that have been found.
#' \item bhattacharrya distance: Bhattacharyya distance between 1-somy and 2-somy distributions.
#' \item sos: Sum-of-squares distance of read counts to the fitted distributions in their respective segments.
#' }
#'
#' @param models A list of \code{\link{GRanges}} or \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @return A data.frame with columns
#' @author Aaron Taudt
#' @export
#' @examples 
#'## Get a list of HMMs
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'df <- getQC(files)
getQC <- function(models) {

    ## Helper function
    null2na <- function(x) {
        if (is.null(x)) {
            return(NA)
        } else {
            return(x)
        }
    }
  	models <- suppressMessages( loadFromFiles(models, check.class=c('GRanges', class.univariate.hmm, class.bivariate.hmm)) )
  	qframe <- list()
  	for (i1 in 1:length(models)) {
    		model <- models[[i1]]
    		if (class(model) == class.univariate.hmm | class(model) == class.bivariate.hmm) {
    		    bins <- model$bins
        		qframe[[i1]] <- data.frame( total.read.count=sum(bins$counts),
                                				avg.binsize=mean(width(bins)),
                                				avg.read.count=mean(bins$counts),
                                				spikiness=qc.spikiness(bins$counts),
                                				entropy=qc.entropy(bins$counts),
                                				complexity=null2na(model$qualityInfo$complexity[1]),
                                				loglik=null2na(model$convergenceInfo$loglik),
                                				num.segments=length(model$segments),
                                				bhattacharyya=qc.bhattacharyya(model),
                                				sos=qc.sos(model)
                                				)
    		} else if (class(model) == 'GRanges') {
    		    bins <- model
        		qframe[[i1]] <- data.frame( total.read.count=sum(bins$counts),
                                				avg.binsize=mean(width(bins)),
                                				avg.read.count=mean(bins$counts),
                                				spikiness=qc.spikiness(bins$counts),
                                				entropy=qc.entropy(bins$counts),
                                				complexity=null2na(attr(bins,'qualityInfo')$complexity[1]),
                                				loglik=NA,
                                				num.segments=NA,
                                				bhattacharyya=NA,
                                				sos=NA
                                				)
    		}
  	}
  	names(qframe) <- names(models)
  	qframe <- do.call(rbind, qframe)
  	return(qframe)

}
												
#' Cluster based on quality variables
#'
#' This function uses the \pkg{\link{mclust}} package to cluster the input samples based on various quality measures.
#'
#' Please see \code{\link{getQC}} for a brief description of the quality measures.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param G An integer vector specifying the number of clusters that are compared. See \code{\link[mclust:Mclust]{Mclust}} for details.
#' @param itmax The maximum number of outer and inner iterations for the \code{\link[mclust:Mclust]{Mclust}} function. See \code{\link[mclust:emControl]{emControl}} for details.
#' @param measures The quality measures that are used for the clustering. Supported is any combination of \code{c('spikiness','entropy','num.segments','bhattacharyya','loglik','complexity','sos','avg.read.count','total.read.count','avg.binsize')}. 
#' @param orderBy The quality measure to order the clusters by. Default is \code{'spikiness'}.
#' @param reverseOrder Logical indicating whether the ordering by \code{orderBy} is reversed.
#' @return A \code{list} with the classification, parameters and the \code{\link[mclust]{Mclust}} fit.
#' @author Aaron Taudt
#' @seealso \code{\link{getQC}}
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
clusterByQuality <- function(hmms, G=1:9, itmax=c(100,100), measures=c('spikiness','entropy','num.segments','bhattacharyya','complexity','sos'), orderBy='spikiness', reverseOrder=FALSE) {
	
	hmms <- loadFromFiles(hmms, check.class=c(class.univariate.hmm, class.bivariate.hmm))
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
