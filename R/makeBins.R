

#' Make fixed-width bins
#'
#' Make fixed-width bins based on given bin size.
#'
#' @param bamfile A BAM file from which the header is read to determine the chromosome lengths. If a \code{bamfile} is specified, option \code{assembly} is ignored.
#' @param assembly An assembly from which the chromosome lengths are determined. Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. This option is ignored if \code{bamfile} is specified.
#' @param chrom.lengths A named character vector with chromosome lengths. Names correspond to chromosomes.
#' @param chromosome.format A character specifying the format of the chromosomes if \code{assembly} is specified. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...). If a \code{bamfile} is supplied, the format will be chosen automatically.
#' @param binsizes A vector of bin sizes in base pairs.
#' @param chromosomes A subset of chromosomes for which the bins are generated.
#' @return A \code{list()} of \code{\link{GRanges}} objects with fixed-width bins.
#' @author Aaron Taudt
#' @importFrom Rsamtools scanBamHeader
#' @export
fixedWidthBins <- function(bamfile=NULL, assembly=NULL, chrom.lengths=NULL, chromosome.format='NCBI', binsizes=1e6, chromosomes=NULL) {

	### Check user input ###
	if (is.null(bamfile) & is.null(assembly) & is.null(chrom.lengths)) {
		stop("Please specify either a 'bamfile', 'assembly' or 'chrom.lengths'")
	}

	### Get chromosome lengths ###
	if (!is.null(bamfile)) {
		ptm <- startTimedMessage(paste0("Reading header from ", bamfile, " ..."))
		file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
		chrom.lengths <- file.header$targets
		stopTimedMessage(ptm)
	} else if (!is.null(assembly)) {
		ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
		df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)
		chrom.lengths <- df$UCSC_seqlengths
		if (chromosome.format=='UCSC') {
			names(chrom.lengths) <- df$UCSC_seqlevels
		} else if (chromosome.format=='NCBI') {
			names(chrom.lengths) <- df$NCBI_seqlevels
		}
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		stopTimedMessage(ptm)
	} else if (!is.null(chrom.lengths)) {
	}
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring)
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}

	### Making fixed-width bins ###
	bins.list <- list()
	for (binsize in binsizes) {
		ptm <- startTimedMessage("Making fixed-width bins for bin size ", binsize, " ...")
		bins <- GenomicRanges::GRangesList()
		skipped.chroms <- character()
		## Loop over chromosomes
		for (chromosome in chroms2use) {
			## Check last incomplete bin
			incomplete.bin <- chrom.lengths[chromosome] %% binsize > 0
			if (incomplete.bin) {
				numbin <- floor(chrom.lengths[chromosome]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
			} else {
				numbin <- chrom.lengths[chromosome]/binsize
			}
			if (numbin == 0) {
				skipped.chroms[chromosome] <- chromosome
				next
			}
			## Initialize vectors
			chroms <- rep(chromosome,numbin)
			reads <- rep(0,numbin)
			start <- seq(from=1, by=binsize, length.out=numbin)
			end <- seq(from=binsize, by=binsize, length.out=numbin)
	# 		end[length(end)] <- seqlengths(data)[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired

			## Create binned chromosome as GRanges object
			bins.chr <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbin),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbin)
							)
			suppressWarnings(
				bins[[chromosome]] <- bins.chr
			)

		}
		## end loop chromosomes

		### Concatenate all chromosomes
		bins <- unlist(bins)
		names(bins) <- NULL
		seqlengths(bins) <- as.integer(chrom.lengths[names(seqlengths(bins))])
		bins.list[[as.character(binsize)]] <- bins
		stopTimedMessage(ptm)

		if (length(skipped.chroms)>0) {
			warning("The following chromosomes were skipped because they are smaller than the binsize: ", paste0(skipped.chroms, collapse=', '))
		}

	}

	return(bins.list)

}


#' Make variably sized bins
#' 
#' Make variable-width bins based on a reference BAM file. This can be a simulated file (produced by \code{\link{simulateReads}} and aligned with your favourite aligner) or a real reference.
#' 
#' Variable-width bins are produced by first binning the reference BAM file with fixed-width bins and selecting the desired number of reads per bin as the (non-zero) maximum of the histogram. A new set of bins is then generated such that every bin contains the desired number of reads.
#' 
#' @param reads A \code{\link{GRanges}} with reads. See \code{\link{bam2GRanges}} and \code{\link{bed2GRanges}}.
#' @param binsizes A vector with binsizes. Resulting bins will be close to the specified binsizes.
#' @return A \code{list()} of \code{\link{GRanges}} objects with variable-width bins.
#' @author Aaron Taudt
#' @export
variableWidthBins <- function(reads, binsizes, chromosomes=NULL) {
	
	### Check user input ###
	chroms.in.data <- seqlevels(reads)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring)
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}

	## Drop unwanted seqlevels
	reads <- reads[seqnames(reads) %in% chroms2use]
	reads <- keepSeqlevels(reads, chroms2use)

	## Make fixed width bins
	ptm <- startTimedMessage(paste0("Binning ", bamfile, " ..."))
	binned.list <- suppressMessages( align2binned(reads, format='GRanges', binsizes=binsizes, calc.complexity=FALSE, chromosomes=chromosomes) )
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


#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param bamfile A sorted BAM file.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @param what A character vector of fields that are returned. Type \code{\link[Rsamtools]{scanBamWhat}} to see what is available.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairsFromBam readGAlignmentsFromBam first
bam2GRanges <- function(bamfile, bamindex=bamfile, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000, what='mapq') {

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		bamindex.own <- Rsamtools::indexBam(bamfile)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
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
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, flag=scanBamFlag(isDuplicate=FALSE)))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, flag=scanBamFlag(isDuplicate=FALSE)))
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

	return(data)

}


#' Import BED file into GRanges
#'
#' Import aligned reads from a BED file into a \code{\link{GRanges}} object.
#'
#' @param bedfile A file with aligned reads in BED format. The columns have to be c('chromosome','start','end','description','mapq','strand').
#' @param assembly Please see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments.
#' @importFrom GenomeInfoDb fetchExtendedChromInfoFromUCSC
bed2GRanges <- function(bedfile, assembly, chromosomes=NULL, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000) {

	# File with reads, specify classes for faster import (0-based)
	ptm <- startTimedMessage("Reading file ",basename(bedfile)," ...")
	classes <- c('character','numeric','numeric','NULL','integer','character')
	data.raw <- read.table(bedfile, colClasses=classes)
	# Convert to GRanges object
	data <- GenomicRanges::GRanges(seqnames=Rle(data.raw[,1]), ranges=IRanges(start=data.raw[,2]+1, end=data.raw[,3]), strand=data.raw[,5])	# start+1 to go from [0,x) -> [1,x]
	mcols(data)$mapq <- data.raw[,4]
	remove(data.raw)
	chroms.in.data <- seqlevels(data)
	stopTimedMessage(ptm)
	# Get chromosome lengths
	ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
	df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)
	chrom.lengths <- df$UCSC_seqlengths
	if (grepl('^chr',seqlevels(data)[1])) {
		names(chrom.lengths) <- df$UCSC_seqlevels
	} else {
		names(chrom.lengths) <- df$NCBI_seqlevels
	}
	seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	stopTimedMessage(ptm)

	## Issue warnings and stuff for non-existing chromosomes
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
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

	return(data)

}


