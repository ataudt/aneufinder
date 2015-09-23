

#' Read bed-file into GRanges
#'
#' This is a simple convenience function to read a (compressed) bed-file into a \code{\link{GRanges}} object. The bed-file is expected to have the following fields: \code{chrom, chromStart, chromEnd, name, score, strand}.
#'
#' @param bedfile Filename of the bed or bed.gz file.
#' @param skip Number of lines to skip at the beginning.
#' @return A \code{\link{GRanges}} object with the contents of the bed-file.
#' @author Aaron Taudt
#' @export
bed2GRanges <- function(bedfile, skip=0) {

	# File with reads, determine classes first for faster import (0-based)
	tab5rows <- read.table(bedfile, nrows=5, skip=skip)
	classes.in.bed <- sapply(tab5rows, class)
	data <- read.table(bedfile, colClasses=classes.in.bed, skip=skip)
	# GRanges compatible strand information
	data[,6] <- sub('.','*',data[,6])
	# Convert to GRanges object
	gr <- GenomicRanges::GRanges(seqnames=Rle(data[,1]),
																	ranges=IRanges(start=data[,2]+1, end=data[,3]+1),	# +1 to match coordinate systems
																	strand=data[,6],
																	name=data[,4],
																	score=data[,5])

	return(gr)

}

