

# =================================================================
# Extraction of segments and clustering
# =================================================================
#' Extract segments and cluster
#'
#' Extract segments and ID from a list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects and cluster if desired.
#'
#' @param hmms A list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects or a character vector of files that contains such objects.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param classes A vector with class labels the same length as \code{hmms}. If supplied, the clustering will be ordered optimally with respect to the class labels (see \code{\link[ReorderCluster]{RearrangeJoseph}}).
#' @param exclude.regions A \code{\link{GRanges}} with regions that will be excluded from the computation of the clustering. This can be useful to exclude regions with artifacts.
#' @return A \code{list()} with (clustered) segments and breakpoint coordinates.
#' @importFrom ReorderCluster RearrangeJoseph
#' @importFrom stats as.dist cov.wt hclust
getSegments <- function(hmms, cluster=TRUE, classes=NULL, exclude.regions=NULL) {

	## Load the files
	hmms <- loadFromFiles(hmms, check.class=c(class.univariate.hmm, class.bivariate.hmm))

	## Get segments from list
	ptm <- startTimedMessage("Getting segments ...")
	segs <- GRangesList()
	hmms2use <- numeric()
	for (i1 in 1:length(hmms)) {
	  hmm <- hmms[[i1]]
		if (!is.null(hmm$segments)) {
		  hmms2use[hmm$ID] <- i1
			segs[[as.character(hmm$ID)]] <- hmm$segments
		}
	}
	hmms <- hmms[hmms2use]
	stopTimedMessage(ptm)

	## Clustering based on bins
	if (cluster) {
		ptm <- startTimedMessage("Making consensus template ...")
		if (!is.null(hmms[[1]]$bins$copy.number)) {
  		constates <- sapply(hmms, function(hmm) { hmm$bins$copy.number })
		} else if (!is.null(hmms[[1]]$bins$state)) {
		  constates <- sapply(hmms, function(hmm) { suppressWarnings( initializeStates(levels(hmm$bins$state))$multiplicity[as.character(hmm$bins$state)] ) })
		}
		constates[is.na(constates)] <- 0
		vars <- apply(constates, 1, var, na.rm=TRUE)
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Clustering ...")
    ## Exclude regions ##
    if (!is.null(exclude.regions)) {
        ind <- findOverlaps(hmms[[1]]$bins, exclude.regions)@from
    		constates <- constates[-ind,]
    }
		dist <- stats::dist(t(constates))
		hc <- stats::hclust(dist)
		stopTimedMessage(ptm)
		# Dendrogram
		message("Reordering ...")
		if (!is.null(classes)) {
			# Reorder by classes
			res <- ReorderCluster::RearrangeJoseph(hc, as.matrix(dist), class=classes, cpp=TRUE)
			file.remove('A.txt','minI.txt','minJ.txt')
			hc <- res$hcl
		}
		# Reorder samples
		segs <- segs[hc$order]

		return(list(segments=segs, clustering=hc, dist=dist))
	}

	return(list(segments=segs))
}
