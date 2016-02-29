

#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param bamfile A sorted BAM file.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your BAM files (not implemented for BED files).
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @param blacklist A \code{\link{GRanges}} or a bed(.gz) file with blacklisted regions. Reads falling into those regions will be discarded.
#' @param what A character vector of fields that are returned. Type \code{\link[Rsamtools]{scanBamWhat}} to see what is available.
#' @return A \code{\link{GRanges}} object containing the reads.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first
#' @export
#'
#'@examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata/example.bam", package="aneufinder")
#'## Read the file into a GRanges object
#'reads <- bam2GRanges(bamfile, chromosomes=c(1:19,'X','Y'), pairedEndReads=FALSE,
#'                     min.mapq=10, remove.duplicate.reads=TRUE)
#'print(reads)
#'
bam2GRanges <- function(bamfile, bamindex=bamfile, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, blacklist=NULL, what='mapq') {

	## Input checks
	if (!is.null(blacklist)) {
		if ( !(is.character(blacklist) | class(blacklist)=='GRanges') ) {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
	}

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		ptm <- startTimedMessage("Making bam-index file ...")
		bamindex.own <- Rsamtools::indexBam(bamfile)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
		stopTimedMessage(ptm)
	}
	file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
	chrom.lengths <- file.header$targets
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data.')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}

	## Import the file into GRanges
	ptm <- startTimedMessage("Reading file ",basename(bamfile)," ...")
	gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, flag=scanBamFlag(isDuplicate=FALSE)))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, flag=scanBamFlag(isDuplicate=FALSE)))
		}
	}
	stopTimedMessage(ptm)

	if (length(data.raw) == 0) {
		if (pairedEndReads) {
			stop(paste0("No reads imported. Does your file really contain paired end reads? Try with 'pairedEndReads=FALSE'"))
		}
		stop(paste0('No reads imported! Check your BAM-file ', bamfile))
	}

	## Filter by mapping quality
	if (pairedEndReads) {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- as(data.raw, 'GRanges') # treat as one fragment
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		if (!is.null(min.mapq)) {
			mapq.first <- mcols(GenomicAlignments::first(data.raw))$mapq
			mapq.last <- mcols(GenomicAlignments::last(data.raw))$mapq
			mapq.mask <- mapq.first >= min.mapq & mapq.last >= min.mapq
			if (any(is.na(mapq.mask))) {
				warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
			}
			data <- data[which(mapq.mask)]
		}
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	} else {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- as(data.raw, 'GRanges')
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		if (!is.null(min.mapq)) {
			if (any(is.na(mcols(data)$mapq))) {
				warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
				mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
			}
			data <- data[mcols(data)$mapq >= min.mapq]
		}
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	}

	## Exclude reads falling into blacklisted regions
	if (!is.null(blacklist)) {
		ptm <- startTimedMessage("Filtering blacklisted regions ...")
		if (is.character(blacklist)) {
			if (grepl('^chr', seqlevels(data)[1])) {
				chromosome.format <- 'UCSC'
			} else {
				chromosome.format <- 'NCBI'
			}
			black <- importBed(blacklist, skip=0, chromosome.format=chromosome.format)
		} else if (class(blacklist)=='GRanges') {
			black <- blacklist
		} else {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
		overlaps <- findOverlaps(data, black)
		idx <- setdiff(1:length(data), queryHits(overlaps))
		data <- data[idx]
		stopTimedMessage(ptm)
	}

	return(data)

}


#' Import BED file into GRanges
#'
#' Import aligned reads from a BED file into a \code{\link{GRanges}} object.
#'
#' @param bedfile A file with aligned reads in BED format. The columns have to be c('chromosome','start','end','description','mapq','strand').
#' @param assembly Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. Only necessary when importing BED files. BAM files are handled automatically. Alternatively a data.frame generated by \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}}.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments.
#' @param blacklist A \code{\link{GRanges}} or a bed(.gz) file with blacklisted regions. Reads falling into those regions will be discarded.
#' @return A \code{\link{GRanges}} object containing the reads.
#' @importFrom utils read.table
#' @export
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata/KK150311-VI_07.bam.bed.gz", package="aneufinder")
#'## Read the file into a GRanges object
#'reads <- bed2GRanges(bedfile, assembly='mm10', chromosomes=c(1:19,'X','Y'),
#'                     min.mapq=10, remove.duplicate.reads=TRUE)
#'print(reads)
#'
bed2GRanges <- function(bedfile, assembly, chromosomes=NULL, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, blacklist=NULL) {

	## Input checks
	if (!is.null(blacklist)) {
		if ( !(is.character(blacklist) | class(blacklist)=='GRanges') ) {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
	}

	# File with reads, specify classes for faster import (0-based)
	ptm <- startTimedMessage("Reading file ",basename(bedfile)," ...")
	classes <- c('character','numeric','numeric','NULL','integer','character')
	data.raw <- utils::read.table(bedfile, colClasses=classes)
	# Convert to GRanges object
	data <- GenomicRanges::GRanges(seqnames=data.raw[,1], ranges=IRanges(start=data.raw[,2]+1, end=data.raw[,3]), strand=data.raw[,5])	# start+1 to go from [0,x) -> [1,x]
	mcols(data)$mapq <- data.raw[,4]
	remove(data.raw)
	stopTimedMessage(ptm)
	# Get chromosome lengths
	if (is.character(assembly)) {
		ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
		df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)
		stopTimedMessage(ptm)
	} else if (is.data.frame(assembly)) {
		df <- assembly
	} else {
		stop("Unknown assembly")
	}
	chrom.lengths <- df$UCSC_seqlength
	if (grepl('^chr',seqlevels(data)[1])) {
		names(chrom.lengths) <- df$UCSC_seqlevel
	} else {
		names(chrom.lengths) <- df$NCBI_seqlevel
	}
	seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])

	chroms.in.data <- seqlevels(data)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if none of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data.')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}

	## Select only desired chromosomes
	data <- data[seqnames(data) %in% chroms2use]
	data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
	## Drop seqlevels where seqlength is NA
	na.seqlevels <- seqlevels(data)[is.na(seqlengths(data))]
	data <- data[seqnames(data) %in% seqlevels(data)[!is.na(seqlengths(data))]]
	data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
	if (length(na.seqlevels) > 0) {
		warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
	}

	if (length(data) == 0) {
		stop(paste0('No reads imported!'))
	}

	## Filter by mapping quality
	ptm <- startTimedMessage("Filtering reads ...")
	if (!is.null(min.mapq)) {
		if (any(is.na(mcols(data)$mapq))) {
			warning(paste0(bedfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
			mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
		}
		data <- data[mcols(data)$mapq >= min.mapq]
	}
	# Filter out too long fragments
	data <- data[width(data)<=max.fragment.width]
	stopTimedMessage(ptm)

	## Exclude reads falling into blacklisted regions
	if (!is.null(blacklist)) {
		ptm <- startTimedMessage("Filtering blacklisted regions ...")
		if (is.character(blacklist)) {
			if (grepl('^chr', seqlevels(data)[1])) {
				chromosome.format <- 'UCSC'
			} else {
				chromosome.format <- 'NCBI'
			}
			black <- importBed(blacklist, skip=0, chromosome.format=chromosome.format)
		} else if (class(blacklist)=='GRanges') {
			black <- blacklist
		} else {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
		overlaps <- findOverlaps(data, black)
		idx <- setdiff(1:length(data), queryHits(overlaps))
		data <- data[idx]
		stopTimedMessage(ptm)
	}

	return(data)

}

