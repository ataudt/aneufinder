

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
#' @importFrom stats dnbinom dbinom dpois
qc.bhattacharyya <- function(hmm) {
	if (class(hmm)=='aneuHMM') {
		distr <- hmm$distributions
	} else if (class(hmm)=='aneuBiHMM') {
		distr <- hmm$distributions$minus
	}
  if (is.null(distr)) {
    return(NA)
  }
  distr <- distr[!rownames(distr) %in% c('zero-inflation', '0-somy'),]
  if (nrow(distr) <= 1) {
    warning(hmm$ID, ": Couldn't calculate Bhattacharyya distance.")
    return(NA)
  }
  if (!is.null(hmm$bins$counts)) {
    counts <- hmm$bins$counts
  } else if (!is.null(hmm$bincounts[[1]]$counts)) {
    counts <- hmm$bincounts[[1]]$counts
  }
	x <- 0:max(max(counts),500)
	if (!is.na(distr['1-somy','size']) & !is.na(distr['2-somy','size'])) {
	  if (distr['1-somy','type'] == 'dnbinom') {
  	  term1 <- stats::dnbinom(x, size=distr['1-somy','size'], prob=distr['1-somy','prob'])
	  } else if (distr['1-somy','type'] == 'dbinom') {
  	  term1 <- stats::dbinom(x, size=round(distr['1-somy','size']), prob=distr['1-somy','prob'])
	  } else if (distr['1-somy','type'] == 'dpois') {
  	  term1 <- stats::dpois(x, lambda=distr['1-somy','prob'])
	  }
	  if (distr['2-somy','type'] == 'dnbinom') {
  	  term2 <- stats::dnbinom(x, size=distr['2-somy','size'], prob=distr['2-somy','prob'])
	  } else if (distr['2-somy','type'] == 'dbinom') {
  	  term2 <- stats::dbinom(x, size=round(distr['2-somy','size']), prob=distr['2-somy','prob'])
	  } else if (distr['2-somy','type'] == 'dpois') {
  	  term2 <- stats::dpois(x, lambda=distr['2-somy','prob'])
	  }
	} else {
	  if (distr[1,'type'] == 'dnbinom') {
  	  term1 <- stats::dnbinom(x, size=distr[1,'size'], prob=distr[1,'prob'])
	  } else if (distr[1,'type'] == 'dbinom') {
  	  term1 <- stats::dbinom(x, size=round(distr[1,'size']), prob=distr[1,'prob'])
	  } else if (distr[1,'type'] == 'dpois') {
  	  term1 <- stats::dpois(x, lambda=distr[1,'prob'])
	  }
	  if (distr[2,'type'] == 'dnbinom') {
  	  term2 <- stats::dnbinom(x, size=distr[2,'size'], prob=distr[2,'prob'])
	  } else if (distr[2,'type'] == 'dbinom') {
  	  term2 <- stats::dbinom(x, size=round(distr[2,'size']), prob=distr[2,'prob'])
	  } else if (distr[2,'type'] == 'dpois') {
  	  term2 <- stats::dpois(x, lambda=distr[2,'prob'])
	  }
  	warning(hmm$ID, ": Bhattacharyya distance calculated for ", rownames(distr)[1], " and ", rownames(distr)[2], " instead of 1-somy and 2-somy.")
	}
	dist <- -log(sum(sqrt(term1 * term2)))
	return(dist)
}

#' @describeIn qualityControl Sum-of-squares distance from the read counts to the fitted distributions
qc.sos <- function(hmm) {
  if (!is.null(hmm$bins$counts)) {
    counts <- hmm$bins$counts
    mcounts <- hmm$bins$mcounts
    pcounts <- hmm$bins$pcounts
    state <- hmm$bins$state
    mstate <- hmm$bins$mstate
    pstate <- hmm$bins$pstate
  } else if (!is.null(hmm$bincounts[[1]]$counts)) {
    counts <- hmm$bincounts[[1]]$counts
    mcounts <- hmm$bincounts[[1]]$mcounts
    pcounts <- hmm$bincounts[[1]]$pcounts
    ind <- findOverlaps(hmm$bins, hmm$bincounts[[1]], select='first')
    state <- hmm$bins$state[ind]
    mstate <- hmm$bins$mstate[ind]
    pstate <- hmm$bins$pstate[ind]
  }
	if (class(hmm)=='aneuHMM') {
	  state <- hmm$bins$state
		distr <- hmm$distributions
    mu <- distr$mu
    names(mu) <- rownames(distr)
    sos <- sum( (counts - mu[as.character(state)]) ^ 2 )
	} else if (class(hmm)=='aneuBiHMM') {
	  state <- 
		distr <- hmm$distributions$minus
    mu <- distr$mu
    names(mu) <- rownames(distr)
    sos <- sum( (c(mcounts - mu[as.character(mstate)],
                   pcounts - mu[as.character(pstate)])
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
#' @param models A list of \code{\link{GRanges-class}} or \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
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
  	models <- suppressMessages( loadFromFiles(models, check.class=c('GRanges', 'GRangesList', "aneuHMM", "aneuBiHMM")) )
  	qframe <- list()
  	for (i1 in 1:length(models)) {
    		model <- models[[i1]]
    		if (class(model) == "aneuHMM" | class(model) == "aneuBiHMM") {
    		    bins <- model$bins
            if (!is.null(model$bins$counts)) {
              counts <- model$bins$counts
            } else if (!is.null(model$bincounts[[1]]$counts)) {
              counts <- model$bincounts[[1]]$counts
            }
        		qframe[[i1]] <- data.frame( total.read.count=sum(counts),
                                				avg.binsize=mean(width(bins)),
                                				avg.read.count=mean(counts),
                                				spikiness=qc.spikiness(counts),
                                				entropy=qc.entropy(counts),
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
    		} else if (is(model, "GRangesList")) {
    		    bins <- model[[1]]
        		qframe[[i1]] <- data.frame( total.read.count=sum(bins$counts),
                                				avg.binsize=mean(width(bins)),
                                				avg.read.count=mean(bins$counts),
                                				spikiness=qc.spikiness(bins$counts),
                                				entropy=qc.entropy(bins$counts),
                                				complexity=null2na(attr(model,'qualityInfo')$complexity[1]),
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
	
	hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
	df <- getQC(hmms)
	df <- df[measures]
	ptm <- startTimedMessage("clustering ...")
	# Remove NA/NaN/Inf columns
	df[apply(df, 2, is.infinite)] <- NA
	df[apply(df, 2, is.nan)] <- NA
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
