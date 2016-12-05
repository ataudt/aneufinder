

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
#' @return A \code{list()} with (clustered) segments and SCE coordinates.
#' @importFrom ReorderCluster RearrangeJoseph
#' @importFrom stats as.dist cov.wt hclust
getSegments <- function(hmms, cluster=TRUE, classes=NULL) {

	## Load the files
	hmms <- loadFromFiles(hmms, check.class=c(class.univariate.hmm, class.bivariate.hmm))

	## Get segments from list
	ptm <- startTimedMessage("Getting segments ...")
	grlred <- GRangesList()
	hmms2use <- numeric()
	for (i1 in 1:length(hmms)) {
	  hmm <- hmms[[i1]]
		if (!is.null(hmm$segments)) {
		  hmms2use[hmm$ID] <- i1
			grlred[[as.character(hmm$ID)]] <- hmm$segments
		}
	}
	hmms <- hmms[hmms2use]
	stopTimedMessage(ptm)

	## Clustering based on bins
	if (cluster) {
		ptm <- startTimedMessage("Making consensus template ...")
		constates <- sapply(hmms, function(hmm) { hmm$bins$copy.number })
		constates[is.na(constates)] <- 0
		meanstates <- apply(constates, 1, mean, na.rm=TRUE)
		vars <- apply(constates, 1, var, na.rm=TRUE)
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Clustering ...")
		# Null bins with artificially high variance
		constates[vars >= quantile(vars, 0.999)] <- 0
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
		grlred <- grlred[hc$order]

		return(list(segments=grlred, clustering=hc, dist=dist))
	}

	return(list(segments=grlred))
}
