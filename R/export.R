

#' Export genome browser viewable files
#'
#' Export copy-number-variation state or read counts as genome browser viewable file
#'
#' Use \code{exportCNVs} to export the copy-number-variation state from an \code{\link{aneuHMM}} object in BED format.
#' Use \code{exportReadCounts} to export the binned read counts from an \code{\link{aneuHMM}} object in WIGGLE format.
#' Use \code{exportGRanges} to export a \code{\link{GRanges}} object in BED format.
#'
#' @return \code{NULL}
#' @name export
#' @author Aaron Taudt
#'@examples
#'\dontrun{
#'## Get results from a small-cell-lung-cancer
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'## Export the CNV states for upload to the UCSC genome browser
#'exportCNVs(files, filename='upload-me-to-a-genome-browser', cluster=TRUE)}
#'
NULL


# ==================================
# Replace '1' by 'chr1' if necessary
# ==================================
insertchr <- function(hmm.gr) {
	mask <- which(!grepl('chr', seqnames(hmm.gr)))
	mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
	mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
	mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
	return(hmm.gr)
}
stripchr <- function(hmm.gr) {
	# Change chromosome names from 'chr1' to '1' if necessary
	mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
	mcols(hmm.gr)$chromosome <- sub(pattern='^chr', replacement='', mcols(hmm.gr)$chromosome)
	mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
	return(hmm.gr)
}

# ==============================================
# Write color-coded tracks with states from HMMs
# ==============================================
#' @describeIn export Export CNV-state as .bed.gz file
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for CNV-state or ".wig.gz" for read counts. Any existing file will be overwritten.
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity in their CNV-state.
#' @param export.CNV A logical, indicating whether the CNV-state shall be exported.
#' @param export.SCE A logical, indicating whether breakpoints shall be exported.
#' @importFrom grDevices col2rgb
#' @importFrom utils write.table
#' @export
exportCNVs <- function(hmms, filename, cluster=TRUE, export.CNV=TRUE, export.SCE=TRUE) {

	## Get segments and breakpoint coordinates
	hmms <- loadFromFiles(hmms, check.class=c(class.univariate.hmm, class.bivariate.hmm))
	temp <- getSegments(hmms, cluster=cluster)
	hmm.grl <- temp$segments
	if (cluster) {
		hmms <- hmms[temp$clustering$order]
	}
	if (export.SCE) {
		breakpoints <- lapply(hmms,'[[','breakpoints')
		names(breakpoints) <- lapply(hmms,'[[','ID')
		breakpoints <- breakpoints[!unlist(lapply(breakpoints, is.null))]
		breakpoints <- breakpoints[lapply(breakpoints, length)!=0]		
		if (length(breakpoints)==0) {
			export.SCE <- FALSE
		}
	}
	
	### CNV-state ###
	if (export.CNV) {
		# Replace '1' by 'chr1' if necessary
		hmm.grl <- endoapply(hmm.grl, insertchr)
		# Variables
		nummod <- length(hmm.grl)
		filename.bed <- paste0(filename,"_CNV.bed.gz")
		# Generate the colors
		colors <- stateColors(levels(hmm.grl[[1]]$state))
		RGBs <- t(grDevices::col2rgb(colors))
		RGBs <- apply(RGBs,1,paste,collapse=",")
		# Write first line to file
		message('writing to file ',filename.bed)
		filename.gz <- gzfile(filename.bed, 'w')
		cat("", file=filename.gz)
		
		## Write every model to file
		for (imod in 1:nummod) {
			message('writing hmm ',imod,' / ',nummod)
			hmm.gr <- hmm.grl[[imod]]
			priority <- 51 + 3*imod
			cat(paste0("track name=\"CNV state for ",names(hmm.grl)[imod],"\" description=\"CNV state for ",names(hmm.grl)[imod],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			collapsed.calls <- as.data.frame(hmm.gr)[,c('chromosome','start','end','state')]
			itemRgb <- RGBs[as.character(collapsed.calls$state)]
			numsegments <- nrow(collapsed.calls)
			df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$thickStart - 1
			# Write to file
			utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		}
		close(filename.gz)
	}

	### Breakpoints ###
	if (export.SCE) {
		# Replace '1' by 'chr1' if necessary
		breakpoints <- endoapply(breakpoints, insertchr)
		# Variables
		nummod <- length(breakpoints)
		filename.bed <- paste0(filename,"_SCE.bed.gz")
		# Write first line to file
		message('writing breakpoints to file ',filename.bed)
		filename.gz <- gzfile(filename.bed, 'w')
		cat("", file=filename.gz)
		
		## Write every model to file
		for (imod in 1:nummod) {
			message('writing hmm ',imod,' / ',nummod)
			hmm.gr <- breakpoints[[imod]]
			priority <- 52 + 3*imod
			cat(paste0("track name=\"breakpoints for ",names(breakpoints)[imod],"\" description=\"breakpoints for ",names(breakpoints)[imod],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			collapsed.calls <- as.data.frame(hmm.gr)[,c('chromosome','start','end')]
			collapsed.calls$name <- paste0('breakpoint_',1:nrow(collapsed.calls))
			numsegments <- nrow(collapsed.calls)
			df <- cbind(collapsed.calls, score=0, strand=".", thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end)
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$thickStart - 1
			# Write to file
			utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		}
		close(filename.gz)
	}

}


# =============================
# Write signal tracks from HMMs
# =============================
#' @describeIn export Export binned read counts as .wig.gz file
#' @importFrom utils write.table
#' @export
exportReadCounts <- function(hmms, filename) {

	## Load models
	hmms <- loadFromFiles(hmms, check.class=c(class.univariate.hmm, class.bivariate.hmm))

	## Transform to GRanges
	grl <- lapply(hmms, '[[', 'bins')
	hmm.grl <- endoapply(grl, insertchr)

	# Variables
	nummod <- length(hmms)
	filename <- paste0(filename,".wig.gz")
	filename.gz <- gzfile(filename, 'w')

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod)
		hmm <- hmms[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 50 + 3*imod
		binsize <- width(hmm.gr[1])
		cat(paste0("track type=wig_0 name=\"read count for ",hmm$ID,"\" description=\"read count for ",hmm$ID,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
		# Write read data
		for (chrom in unique(hmm.gr$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			utils::write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
}


#====================================================
# Export regions from GRanges
#====================================================
#' @describeIn export Export \code{\link{GRanges}} object as BED file.
#' @param gr A \code{\link{GRanges}} object.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param trackname The name that will be used as track name and description in the header.
#' @param score A vector of the same length as \code{gr}, which will be used for the 'score' column in the BED file.
#' @param priority Priority of the track for display in the genome browser.
#' @param append Append to \code{filename}.
#' @param chromosome.format A character specifying the format of the chromosomes if \code{assembly} is specified. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).#' @importFrom utils write.table
#' @export
exportGRanges <- function(gr, filename, header=TRUE, trackname=NULL, score=NULL, priority=NULL, append=FALSE, chromosome.format='UCSC') {

	if (header) {
		if (is.null(trackname)) {
			stop("argument 'trackname' must be specified if 'header=TRUE'")
		}
	}
	## Transform to GRanges
	if (chromosome.format=='UCSC') {
		gr <- insertchr(gr)
	} else if (chromosome.format=='NCBI') {
		gr <- stripchr(gr)
	} else {
		stop("Unknown 'chromosome.format'")
	}

	# Variables
	filename <- paste0(filename,".bed.gz")
	if (append) {
		filename.gz <- gzfile(filename, 'a')
	} else {
		filename.gz <- gzfile(filename, 'w')
	}

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz, append=TRUE)
	if (header) {
		strand.colors <- paste0(apply(col2rgb(strandColors(c('+','-'))), 2, function(x) { paste0(x, collapse=',') }), collapse=' ')
		cat(paste0('track name="',trackname,'" description="',trackname,'" visibility=1 colorByStrand="',strand.colors,'" priority=',priority,'\n'), file=filename.gz, append=TRUE)
	}
	if (length(gr)==0) {
  	close(filename.gz)
  	message('')
		return()
	}

	
	### Write model to file ###
	names(gr) <- NULL #delete rownames otherwise as.data.frame can blow up
	regions <- as.data.frame(gr)[c('start','end','strand')]
	regions$name <- 1:nrow(regions)
	regions <- cbind(regions, as.data.frame(mcols(gr)))
	df <- regions[c('chromosome','start','end','name')]
	if (!is.null(score)) {
		df$score <- score
	} else {
		df$score <- 0
	}
	df$strand <- gsub('\\*','.',regions$strand)
	# Convert from 1-based closed to 0-based half open
	df$start <- df$start - 1
	if (nrow(df) == 0) {
		warning('No regions in input')
	} else {
		utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
	}

	close(filename.gz)
	message('')

}



