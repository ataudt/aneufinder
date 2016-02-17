

#' Make variably sized bins
#' 
#' Make variable-width bins based on a reference BAM file. This can be a simulated file (produced by \code{\link{simulateReads}} and aligned with your favourite aligner) or a real reference.
#' 
#' Variable-width bins are produced by first binning the reference BAM file with fixed-width bins and selecting the desired number of reads per bin as the (non-zero) maximum of the histogram. A new set of bins is then generated such that every bin contains the desired number of reads.
#' 
#' @inheritParams bam2GRanges
#' @param binsizes A vector with binsizes. Resulting bins will be close to the specified binsizes.
#' @return A \code{list} with one \code{\link{GRanges}} object for each element in \code{binsizes}.
#' @author Aaron Taudt
#' @export
variableWidthBins <- function(bamfile, bamindex=bamfile, binsizes, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000) {
	
	## Make fixed width bins
	reads <- bam2GRanges(bamfile, bamindex=bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width)
	ptm <- startTimedMessage(paste0("Binning ", bamfile, " ..."))
	binned.list <- suppressMessages( bam2binned(reads, binsizes=binsizes, calc.complexity=FALSE, chromosomes=chromosomes) )
	stopTimedMessage(ptm)
	
	## Loop over binsizes
	bins.list <- list()
	for (i1 in 1:length(binsizes)) {
		binsize <- binsizes[i1]
		binned <- binned.list[[i1]]
		## Get mode of histogram
		tab <- table(binned$counts)
		modecount <- as.integer(names(which.max(tab[names(tab)!=0])))
		## Pick only every modecount read
		idx <- seq(modecount, length(reads), by=modecount)
		subreads <- reads[idx]
		strand(subreads) <- '*'
		## Make new bins
		bins <- gaps(subreads)
		bins <- bins[strand(bins)=='*']
		bins.list[[as.character(binsize)]] <- bins
	}
	
	return(bins.list)

}
