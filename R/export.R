

#' Export genome browser viewable files
#'
#' Export copy-number-variation state or read counts as genome browser viewable file
#'
#' Use \code{exportCNVs} to export the copy-number-variation state from an \code{\link{aneuHMM}} object in BED format.
#' Use \code{exportReadCounts} to export the binned read counts from an \code{\link{aneuHMM}} object in WIGGLE format.
#' Use \code{exportGRanges} to export a \code{\link{GRanges}} object in BED format.
#'
#' @name export
#' @author Aaron Taudt
#'@examples
#'\dontrun{
#'## Get results from a small-cell-lung-cancer
#'folder <- system.file("extdata/primary-lung/results_univariate", package="aneufinder")
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

# ==============================================
# Write color-coded tracks with states from HMMs
# ==============================================
#' @describeIn export Export CNV-state as .bed.gz file
#' @param hmms A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for CNV-state or ".wig.gz" for read counts. Any existing file will be overwritten.
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity in their CNV-state.
#' @param export.CNV A logical, indicating whether the CNV-state shall be exported.
#' @param export.SCE A logical, indicating whether the SCE events shall be exported.
#' @export
exportCNVs <- function(hmms, filename, cluster=TRUE, export.CNV=TRUE, export.SCE=TRUE) {

	## Get segments and SCE coordinates
	hmms <- loadHmmsFromFiles(hmms)
	temp <- getSegments(hmms, cluster=cluster)
	hmm.grl <- temp$segments
	if (cluster) {
		hmms <- hmms[temp$clustering$order]
	}
	if (export.SCE) {
		sce <- lapply(hmms,'[[','sce')
		names(sce) <- lapply(hmms,'[[','ID')
		sce <- sce[!unlist(lapply(sce, is.null))]
		sce <- sce[lapply(sce, length)!=0]		
		if (length(sce)==0) {
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
		colors <- stateColors()[levels(hmm.grl[[1]]$state)]
		RGBs <- t(col2rgb(colors))
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
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
		close(filename.gz)
	}

	### SCE events ###
	if (export.SCE) {
		# Replace '1' by 'chr1' if necessary
		sce <- endoapply(sce, insertchr)
		# Variables
		nummod <- length(sce)
		filename.bed <- paste0(filename,"_SCE.bed.gz")
		# Write first line to file
		message('writing SCE to file ',filename.bed)
		filename.gz <- gzfile(filename.bed, 'w')
		cat("", file=filename.gz)
		
		## Write every model to file
		for (imod in 1:nummod) {
			message('writing hmm ',imod,' / ',nummod)
			hmm.gr <- sce[[imod]]
			priority <- 52 + 3*imod
			cat(paste0("track name=\"SCE events for ",names(sce)[imod],"\" description=\"SCE events for ",names(sce)[imod],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			collapsed.calls <- as.data.frame(hmm.gr)[,c('chromosome','start','end')]
			collapsed.calls$name <- paste0('SCE_',1:nrow(collapsed.calls))
			numsegments <- nrow(collapsed.calls)
			df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end)
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$thickStart - 1
			# Write to file
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
		close(filename.gz)
	}

}


# =============================
# Write signal tracks from HMMs
# =============================
#' @describeIn export Export binned read counts as .wig.gz file
#' @export
exportReadCounts <- function(hmms, filename) {

	## Load models
	hmms <- loadHmmsFromFiles(hmms)

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
			write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
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
#' @export
exportGRanges <- function(gr, filename, header=TRUE, trackname=NULL, score=NULL, priority=NULL, append=FALSE) {

	if (header) {
		if (is.null(trackname)) {
			stop("argument 'trackname' must be specified if 'header=TRUE'")
		}
	}
	if (length(gr)==0) {
		warning("Supplied GRanges object contains no ranges. Nothing exported.")
		return()
	}

	## Function definitions
	insertchr <- function(hmm.gr) {
		# Change chromosome names from '1' to 'chr1' if necessary
		mask <- which(!grepl('chr', seqnames(hmm.gr)))
		mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
		mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
		mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
		return(hmm.gr)
	}

	## Transform to GRanges
	gr <- insertchr(gr)

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
		cat(paste0('track name="',trackname,'" description="',trackname,'" visibility=1 colorByStrand="',paste(strandColors(), collapse=' '),'" priority=',priority,'\n'), file=filename.gz, append=TRUE)
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
		write.table(format(df, scientific=FALSE), file=filename.gz, append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

	close(filename.gz)
	message('')

}



