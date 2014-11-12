cluster.hmms <- function(hmm.list) {

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)
	
	## Transform to GRanges in reduced representation
ptm <- proc.time()
	temp <- hmmList2GRangesList(hmm.list, reduce=T, numCPU=numCPU, consensus=T)
	grlred <- temp$grl
	consensus <- temp$consensus
print(proc.time()-ptm)
	
ptm <- proc.time()
	## Cluster the samples
	cat("clustering the samples ...")
	# Find states along the consensus template
	library(foreach)
	constates <- foreach (gr = grlred, .packages='GenomicRanges', .combine='cbind') %do% {
		splt <- split(gr, mcols(gr)$state)
		mind <- as.matrix(findOverlaps(consensus, splt))
		col <- matrix(-1, nrow=length(consensus), ncol=1)
		col[mind[,'queryHits'],1] <- mind[,'subjectHits']
		col
	}
	colnames(constates) <- unlist(lapply(hmm.list, '[[', 'ID'))
	# Distance measure
	wcor <- cov.wt(constates, wt=as.numeric(width(consensus)), cor=T)
	dist <- as.dist(1-wcor$cor)
	# Dendrogram
	hc <- hclust(dist)
	cat(" done\n")
print(proc.time()-ptm)

	return(hc)
}
