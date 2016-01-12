

# =================================================================
# Extraction of segments and clustering
# =================================================================
#' Extract segments and cluster
#'
#' Extract segments and ID from a list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects and cluster if desired.
#'
#' @param hmms A list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects or files that contain such objects.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @return A \code{list()} with (clustered) segments and SCE coordinates.
#' @importFrom ReorderCluster RearrangeJoseph
getSegments <- function(hmms, cluster=TRUE, classes=NULL) {

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)

	## Get segments from list
	grlred <- GRangesList()
	for (hmm in hmms) {
		if (!is.null(hmm$segments)) {
			grlred[[hmm$ID]] <- hmm$segments
		}
	}

	## Clustering
	if (cluster) {
		message("Making consensus template ...", appendLF=F); ptm <- proc.time()
		consensus <- disjoin(unlist(grlred))
		constates <- matrix(NA, ncol=length(grlred), nrow=length(consensus))
		for (i1 in 1:length(grlred)) {
			grred <- grlred[[i1]]
			splt <- split(grred, mcols(grred)$state)
			mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
			constates[,i1] <- mind
		}
		meanstates <- apply(constates, 1, mean, na.rm=T)
		mcols(consensus)$meanstate <- meanstates
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		# Distance measure
		# Use covariance instead of correlation to avoid NaNs for which the hclust fails with error
		message("clustering ...", appendLF=F); ptm <- proc.time()
		constates[is.na(constates)] <- 0
		wcor <- cov.wt(constates, wt=as.numeric(width(consensus)))
		dist <- as.dist(max(wcor$cov)-wcor$cov)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		# Dendrogram
		message("reordering ...")
		hc <- hclust(dist)
		if (!is.null(classes)) {
			# Reorder by classes
			res <- ReorderCluster::RearrangeJoseph(hc, as.matrix(dist), class=classes, cpp=TRUE)
			file.remove('A.txt','minl.txt','minJ.txt')
			hc <- res$hcl
		}
		# Reorder samples
		grlred <- grlred[hc$order]

		return(list(segments=grlred, clustering=hc, dist=dist))
	}

	return(list(segments=grlred))
}
