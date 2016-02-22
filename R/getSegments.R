

# =================================================================
# Extraction of segments and clustering
# =================================================================
#' Extract segments and cluster
#'
#' Extract segments and ID from a list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects and cluster if desired.
#'
#' @param hmms A list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects or files that contain such objects.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param classes A vector with class labels the same length as \code{hmms}. If supplied, the clustering will be ordered optimally with respect to the class labels (see \code{\link[ReorderCluster]{RearrangeJoseph}}).
#' @return A \code{list()} with (clustered) segments and SCE coordinates.
#' @importFrom ReorderCluster RearrangeJoseph
#' @importFrom stats as.dist cov.wt hclust
getSegments <- function(hmms, cluster=TRUE, classes=NULL) {

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)

	## Get segments from list
	ptm <- startTimedMessage("Getting segments ...")
	grlred <- GRangesList()
	for (hmm in hmms) {
		if (!is.null(hmm$segments)) {
			grlred[[as.character(hmm$ID)]] <- hmm$segments
		}
	}
	stopTimedMessage(ptm)

	## Clustering
	if (cluster) {
		ptm <- startTimedMessage("Making consensus template ...")
		consensus <- disjoin(unlist(grlred))
		constates <- matrix(NA, ncol=length(grlred), nrow=length(consensus))
		for (i1 in 1:length(grlred)) {
			grred <- grlred[[i1]]
			splt <- split(grred, mcols(grred)$state)
			mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
			constates[,i1] <- mind
		}
		meanstates <- apply(constates, 1, mean, na.rm=TRUE)
		mcols(consensus)$meanstate <- meanstates
		stopTimedMessage(ptm)

		# Distance measure
		# Use covariance instead of correlation to avoid NaNs for which the hclust fails with error
		ptm <- startTimedMessage("clustering ...")
		constates[is.na(constates)] <- 0
		wcor <- stats::cov.wt(constates, wt=as.numeric(width(consensus)))
		dist <- stats::as.dist(max(wcor$cov)-wcor$cov)
		stopTimedMessage(ptm)
		# Dendrogram
		message("reordering ...")
		hc <- stats::hclust(dist)
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
