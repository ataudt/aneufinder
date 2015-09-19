# aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
# Copyright (C) 2015  Aaron Taudt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


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

