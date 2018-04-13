

#' Read bed-file into GRanges
#'
#' This is a simple convenience function to read a bed(.gz)-file into a \code{\link{GRanges-class}} object. The bed-file is expected to have the following fields: \code{chromosome, start, end, name, score, strand}.
#'
#' @param bedfile Filename of the bed or bed.gz file.
#' @param skip Number of lines to skip at the beginning.
#' @param chromosome.format Desired format of the chromosomes. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).
#' @return A \code{\link{GRanges-class}} object with the contents of the bed-file.
#' @author Aaron Taudt
#' @importFrom utils read.table
#' @export
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Import the file and skip the first 10 lines
#'data <- importBed(bedfile, skip=10)
#'
importBed <- function(bedfile, skip=0, chromosome.format='NCBI') {

	# File with reads, determine classes first for faster import (0-based)
	classes <- c('character','numeric','numeric','character','integer','character')
	data <- utils::read.table(bedfile, colClasses=classes, skip=skip)
	# GRanges compatible strand information
	data[,6] <- sub('.','*',data[,6])
	# Adjust chromosome format
	data[,1] <- sub('^chr', '', data[,1])
	if (chromosome.format=='UCSC') {
		data[,1] <- paste0('chr', data[,1])
	}
	# Convert to GRanges object
	gr <- GenomicRanges::GRanges(seqnames=data[,1],
																	ranges=IRanges(start=data[,2]+1, end=data[,3]),	# +1 to match coordinate systems
																	strand=data[,6],
																	name=data[,4],
																	score=data[,5])

	return(gr)

}

