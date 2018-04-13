

#' Export genome browser viewable files
#'
#' Export copy-number-variation state or read counts as genome browser viewable file
#'
#' Use \code{exportCNVs} to export the copy-number-variation state from an \code{\link{aneuHMM}} object in BED format.
#' Use \code{exportReadCounts} to export the binned read counts from an \code{\link{aneuHMM}} object in WIGGLE format.
#' Use \code{exportGRanges} to export a \code{\link{GRanges-class}} object in BED format.
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
# #' @param trackname The name that will be used as track name and description in the header. # already described elsewhere
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity in their CNV-state.
#' @param export.CNV A logical, indicating whether the CNV-state shall be exported.
#' @param export.breakpoints A logical, indicating whether breakpoints shall be exported.
#' @importFrom grDevices col2rgb
#' @importFrom utils write.table
#' @importFrom S4Vectors endoapply
#' @export
exportCNVs <- function(hmms, filename, trackname=NULL, cluster=TRUE, export.CNV=TRUE, export.breakpoints=TRUE) {

  if (length(hmms) == 1 & cluster==TRUE) {
    cluster <- FALSE
    warning("Cannot do clustering because only one object was given.")
  }
	hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
	## Cluster
	cl <- clusterHMMs(hmms, cluster=cluster)
	hmms <- hmms[rev(cl$IDorder)] # reverse order to match heatmapGenomewide
	## Get segments
  segments<- GRangesList()
  for (i1 in 1:length(hmms)) {
      hmm <- hmms[[i1]]
      if (is.null(hmm$segments)) {
          segments[[hmm$ID]] <- GRanges()
      } else {
  	      segments[[hmm$ID]] <- hmm$segments
      }
  }
	## Get breakpoints
	if (export.breakpoints) {
	  breakpoints <- GRangesList()
	  for (i1 in 1:length(hmms)) {
        hmm <- hmms[[i1]]
	      if (is.null(hmm$breakpoints)) {
	          breakpoints[[hmm$ID]] <- GRanges()
	      } else {
    	      breakpoints[[hmm$ID]] <- hmm$breakpoints
	      }
	  }
		if (length(breakpoints)==0) {
			export.breakpoints <- FALSE
		}
	}
	
	### CNV-state ###
	if (export.CNV) {
		# Replace '1' by 'chr1' if necessary
		segments <- endoapply(segments, insertchr)
		# Variables
		nummod <- length(segments)
		filename.bed <- paste0(filename,"_CNV.bed.gz")
		# Generate the colors
		colors <- stateColors(levels(segments[[1]]$state))
		RGBs <- t(grDevices::col2rgb(colors))
		RGBs <- apply(RGBs,1,paste,collapse=",")
		# Write first line to file
		ptm <- startTimedMessage('Writing to file ',filename.bed, ' ...')
		filename.gz <- gzfile(filename.bed, 'w')
		cat("", file=filename.gz)
		
		## Write every model to file
		for (imod in 1:nummod) {
			# message('writing hmm ',imod,' / ',nummod)
			hmm.gr <- segments[[imod]]
			priority <- 51 + 3*imod
			if (!is.null(trackname)) {
  			trackline <- paste0(trackname, ", CNV state for ", names(segments)[imod])
			} else {
			  trackline <- paste0("CNV state for ", names(segments)[imod])
			}
			cat(paste0('track name="', trackline, '" description="', trackline,'" visibility=1 itemRgb=On priority=',priority,'\n'), file=filename.gz, append=TRUE)
			df0 <- as.data.frame(hmm.gr)[,c('chromosome','start','end','state')]
			itemRgb <- RGBs[as.character(df0$state)]
			df <- cbind(df0, score=0, strand=".", thickStart=df0$start, thickEnd=df0$end, itemRgb=itemRgb)
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$thickStart - 1
			# Write to file
			utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		}
		close(filename.gz)
		stopTimedMessage(ptm)
	}

	### Breakpoints ###
	if (export.breakpoints) {
		# Replace '1' by 'chr1' if necessary
		breakpoints <- endoapply(breakpoints, insertchr)
		# Variables
		nummod <- length(breakpoints)
		filename.bed <- paste0(filename,"_breakpoints.bed.gz")
		# Write first line to file
		ptm <- startTimedMessage('Writing breakpoints to file ',filename.bed, ' ...')
		filename.gz <- gzfile(filename.bed, 'w')
		cat("", file=filename.gz)
		
		## Write every model to file
		for (imod in 1:nummod) {
			# message('writing hmm ',imod,' / ',nummod)
			hmm.gr <- breakpoints[[imod]]
			priority <- 52 + 3*imod
			if (!is.null(trackname)) {
  			trackline <- paste0(trackname, ", breakpoints for ", names(segments)[imod])
			} else {
			  trackline <- paste0("breakpoints for ", names(segments)[imod])
			}
			cat(paste0('track name="', trackline, '" description="', trackline,'" visibility=1 itemRgb=On priority=',priority,'\n'), file=filename.gz, append=TRUE)
			if (is.null(hmm.gr$start.conf)) {
			    hmm.gr$start.conf <- start(hmm.gr)
			    hmm.gr$end.conf <- end(hmm.gr)
			}
			df0 <- as.data.frame(hmm.gr)[,c('chromosome','start.conf','end.conf','type','start','end')]
			if (nrow(df0) > 0) {
    			df <- cbind(df0[,c('chromosome','start.conf','end.conf','type')], score=0, strand=".", thickStart=df0$start, thickEnd=df0$end)
    			df$rgb <- apply(col2rgb(breakpointColors()[as.character(df$type)]), 2, function(x) { paste0(x, collapse=',') })
    			# Convert from 1-based closed to 0-based half open
    			df$start.conf <- df$start.conf - 1
    			df$thickStart <- df$thickStart - 1
    			# Write to file
    			utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
			}
		}
		close(filename.gz)
		stopTimedMessage(ptm)
	}

}


# =============================
# Write signal tracks from HMMs
# =============================
#' @describeIn export Export binned read counts as .wig.gz file
#' @importFrom utils write.table
#' @importFrom S4Vectors endoapply
#' @export
exportReadCounts <- function(hmms, filename) {

	## Load models
	hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))

	## Transform to GRanges
	grl <- lapply(hmms, '[[', 'bins')
	hmm.grl <- endoapply(grl, insertchr)

	# Variables
	nummod <- length(hmms)
	filename <- paste0(filename,".wig.gz")
	filename.gz <- gzfile(filename, 'w')

	# Write first line to file
	ptm <- startTimedMessage('Writing to file ',filename, ' ...')
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		# message('writing hmm ',imod,' / ',nummod)
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
	stopTimedMessage(ptm)
}


#====================================================
# Export regions from GRanges
#====================================================
#' @describeIn export Export \code{\link{GRanges-class}} object as BED file.
#' @param gr A \code{\link{GRanges-class}} object.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param trackname The name that will be used as track name and description in the header.
#' @param score A vector of the same length as \code{gr}, which will be used for the 'score' column in the BED file.
#' @param priority Priority of the track for display in the genome browser.
#' @param append Append to \code{filename}.
#' @param chromosome.format A character specifying the format of the chromosomes if \code{assembly} is specified. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).#' @importFrom utils write.table
#' @param thickStart,thickEnd A vector of the same length as \code{gr}, which will be used for the 'thickStart' and 'thickEnd' columns in the BED file.
#' @param as.wiggle A logical indicating whether a variableStep-wiggle file will be exported instead of a BED file. If \code{TRUE}, \code{wiggle.value} must be specified.
#' @param wiggle.val A vector of the same length as \code{gr}, which will be used for the values in the wiggle file.
#' @export
exportGRanges <- function(gr, filename, header=TRUE, trackname=NULL, score=NULL, priority=NULL, append=FALSE, chromosome.format='UCSC', thickStart=NULL, thickEnd=NULL, as.wiggle=FALSE, wiggle.val) {

  ## Check input
	if (header) {
		if (is.null(trackname)) {
			stop("Argument 'trackname' must be specified if 'header=TRUE'")
		}
	}
  if (!is.null(thickStart) | !is.null(thickEnd)) {
    if (is.null(thickStart) | is.null(thickEnd)) {
      stop("Both 'thickStart' and 'thickEnd' must be specified.")
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
	if (as.wiggle) {
  	filename <- paste0(filename,".wig.gz")
	} else {
  	filename <- paste0(filename,".bed.gz")
	}
	if (append) {
		filename.gz <- gzfile(filename, 'a')
	} else {
		filename.gz <- gzfile(filename, 'w')
	}

	# Write first line to file
	ptm <- startTimedMessage('Writing to file ',filename, ' ...')
	cat("", file=filename.gz, append=TRUE)
	if (header) {
		strand.colors <- paste0(apply(col2rgb(strandColors(c('+','-'))), 2, function(x) { paste0(x, collapse=',') }), collapse=' ')
		header.string <- paste0('track name="',trackname,'" description="',trackname,'" visibility=1 colorByStrand="',strand.colors,'" priority=',priority,'\n')
		if (as.wiggle) {
  		header.string <- paste0('track type="wiggle_0" name="',trackname,'" description="',trackname,'" visibility=1 priority=',priority,' autoScale=On alwaysZero=On\n')
		}
		cat(header.string, file=filename.gz, append=TRUE)
	}
	if (length(gr)==0) {
  	close(filename.gz)
	  stopTimedMessage(ptm)
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
	if (!is.null(thickStart)) {
	  df$thickStart <- thickStart
	}
	if (!is.null(thickEnd)) {
	  df$thickEnd <- thickEnd
	}
	# Convert from 1-based closed to 0-based half open
	df$start <- df$start - 1
	if (!is.null(df$thickStart)) {
	  df$thickStart <- df$thickStart - 1
	}
	if (nrow(df) == 0) {
		warning('No regions in input')
	} else {
	  if (as.wiggle) {
	    dfwig <- data.frame(start=df$start, value=wiggle.val)
	    mask.inf <- is.infinite(wiggle.val)
	    dfwig$value[mask.inf] <- max(wiggle.val[!mask.inf])
	    for (chrom in unique(df$chromosome)) {
	      cat(paste0("variableStep chrom=", chrom, "\n"), file=filename.gz, append=TRUE)
    		utils::write.table(format(dfwig[df$chromosome==chrom,], scientific=FALSE, trim=TRUE), file=filename.gz, append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
	    }
	  } else {
  		utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
	  }
	}

	close(filename.gz)
  stopTimedMessage(ptm)

}



