# ==============================================
# Write color-coded tracks with states from HMMs
# ==============================================
export.hmm2bed <- function(hmm.list, only.modified=TRUE, filename="view_me_in_genome_browser") {

	## Function definitions
	insertchr <- function(hmm.gr) {
		# Change chromosome names from '1' to 'chr1' if necessary
		mask <- which(!grepl('chr', seqnames(hmm.gr)))
		mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
		mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
		mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
		return(hmm.gr)
	}

	## Load models
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges
	hmm.grl <- lapply(hmm.list, '[[', 'segments')
	hmm.grl <- lapply(hmm.grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	filename <- paste0(filename,".bed.gz")
	filename.gz <- gzfile(filename, 'w')

	# Generate the colors
	colors <- state.colors[levels(hmm.grl[[1]]$state)]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")

	# Write first line to file
	message('writing to file ',filename)
# 	cat("browser hide all\n", file=filename.gz)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 51 + 3*imod
		cat(paste0("track name=\"univariate calls for ",hmm$ID,"\" description=\"univariate calls for ",hmm$ID,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
		collapsed.calls <- as.data.frame(hmm.gr)[c('chromosome','start','end','state')]
		itemRgb <- RGBs[as.character(collapsed.calls$state)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		df$thickStart <- df$thickStart - 1
		if (only.modified) {
			df <- df[df$state=='modified',]
		}
		if (nrow(df) == 0) {
			warning('hmm ',imod,' does not contain any \'modified\' calls')
		} else {
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}
	close(filename.gz)
	message('')

}


# =============================
# Write signal tracks from HMMs
# =============================
export.hmm2wiggle <- function(hmm.list, filename="view_me_in_genome_browser") {

	## Function definitions
	insertchr <- function(hmm.gr) {
		# Change chromosome names from '1' to 'chr1' if necessary
		mask <- which(!grepl('chr', seqnames(hmm.gr)))
		mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
		mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
		mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
		return(hmm.gr)
	}

	## Load models
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges
	grl <- lapply(hmm.list, '[[', 'bins')
	hmm.grl <- lapply(grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	filename <- paste0(filename,".wiggle.gz")
	filename.gz <- gzfile(filename, 'w')

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 50 + 3*imod
		binsize <- width(hmm.gr[1])
		cat(paste0("track type=wiggle_0 name=\"read count for ",hmm$ID,"\" description=\"read count for ",hmm$ID,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
		# Write read data
		for (chrom in unique(hmm.gr$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$reads, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
}

