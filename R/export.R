# ====================================
# Write color-coded tracks with states
# ====================================
export.hmm2bed <- function(hmm.list, file="view_me_in_genome_browser", numCPU=1) {

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
	temp <- hmmList2GRangesList(hmm.list, reduce=T, numCPU=numCPU, consensus=T)
	consensus.gr <- insertchr(temp$consensus)
	hmm.grl <- lapply(temp$grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	numstates <- hmm.list[[1]]$num.states
	file <- paste(file,"bed", sep=".")

	# Generate the colors
	colors <- state.colors[levels(hmm.grl[[1]]$state)]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")
	cf <- colorRamp(colors, space='rgb', interpolate='linear')
	consensusRGBs <- round(cf((mcols(consensus.gr)$meanstate-1)/(numstates-1)), 0)
	itemconsensusRgb <- apply(consensusRGBs,1,paste,collapse=",")

	# Write first line to file
	cat('writing to file',file,'\n')
# 	cat("browser hide all\n", file=file)
	cat("", file=file)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		cat('writing hmm',imod,'/',nummod,'\r')
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 51 + imod
		cat(paste0("track name=",hmm$ID," state"," description=\"",hmm$ID," state","\" visibility=1 itemRgb=On priority=",priority,"\n"), file=file, append=TRUE)
		collapsed.calls <- as.data.frame(hmm.gr)[c('chromosome','start','end','state')]
		itemRgb <- RGBs[as.character(collapsed.calls$state)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		write.table(format(df, scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	cat('\n')

	### Write consensus model to file ###
	cat('writing mean states to file\n')
	priority <- 50
	cat(paste0("track name=\"","average state","\" description=\"","average state","\" visibility=1 itemRgb=On priority=",priority,"\n"), file=file, append=TRUE)
	collapsed.calls <- as.data.frame(consensus.gr)[c('chromosome','start','end','meanstate')]
	numsegments <- nrow(collapsed.calls)
	df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemconsensusRgb)
	write.table(format(df, scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

}


# =============================
# Write signal tracks from hmms
# =============================
export.hmm2wiggle <- function(hmm.list, file="view_me_in_genome_browser", numCPU=1) {

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
	temp <- hmmList2GRangesList(hmm.list, reduce=F, numCPU=numCPU, consensus=T)
	hmm.grl <- temp$grl
	consensus.gr <- insertchr(temp$consensus)
	hmm.grl <- lapply(hmm.grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	numstates <- hmm.list[[1]]$num.states
	file <- paste(file,"wiggle", sep=".")

	# Write first line to file
	cat('writing to file',file,'\n')
# 	cat("browser hide all\n", file=file)
	cat("", file=file)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		cat('writing hmm',imod,'/',nummod,'\r')
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 50 + imod
		binsize <- width(hmm.gr[1])
		cat(paste("track type=wiggle_0 name=\"",hmm$ID," read count","\" description=\"",hmm$ID," read count","\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n", sep=""), file=file, append=TRUE)
		# Write read data
		for (chrom in unique(hmm.gr$chromosome)) {
			cat(paste("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n", sep=""), file=file, append=TRUE)
			write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$reads, file=file, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	cat('\n')
}


