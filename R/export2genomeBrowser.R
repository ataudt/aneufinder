# ===============================================
# Write color-coded tracks with univariate states
# ===============================================
univariate2bed <- function(uni.hmm.list, file="view_me_in_genome_browser", numCPU=1) {

	## Intercept user input
	if (check.univariate.modellist(uni.hmm.list)!=0) {
		cat("Loading univariate HMMs from files ...")
		mlist <- NULL
		for (modelfile in uni.hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		uni.hmm.list <- mlist
		remove(mlist)
		cat(" done\n")
		if (check.univariate.modellist(uni.hmm.list)!=0) stop("argument 'uni.hmm.list' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}

	## Transform to GRanges
	cat('transforming to GRanges\n')
	if (numCPU > 1) {
		library(doParallel)
		cl <- makeCluster(numCPU)
		registerDoParallel(cl)
		cfun <- function(...) { GRangesList(...) }
		uni.hmm.grl <- foreach (uni.hmm = uni.hmm.list, .packages='aneufinder', .combine='cfun', .multicombine=TRUE) %dopar% {
			hmm2GRanges(uni.hmm, reduce=T)
		}
		consensus.gr <- disjoin(unlist(uni.hmm.grl))
		constates <- foreach (uni.hmm.gr = uni.hmm.grl, .packages='GenomicRanges', .combine='cbind') %dopar% {
			splt <- split(uni.hmm.gr, mcols(uni.hmm.gr)$state)
			mind <- as.matrix(findOverlaps(consensus.gr, splt, select='first'))
		}
		stopCluster(cl)
	} else {
		uni.hmm.grl <- GRangesList()
		for (uni.hmm in uni.hmm.list) {
			uni.hmm.grl[[length(uni.hmm.grl)+1]] <- hmm2GRanges(uni.hmm, reduce=T)
		}
		consensus.gr <- disjoin(unlist(uni.hmm.grl))
		constates <- matrix(NA, ncol=length(uni.hmm.grl), nrow=length(uni.hmm.grl[[1]]))
		for (i1 in 1:length(uni.hmm.grl)) {
			uni.hmm.gr <- uni.hmm.grl[[i1]]
			splt <- split(uni.hmm.gr, mcols(uni.hmm.gr)$state)
			mind <- as.matrix(findOverlaps(consensus.gr, splt, select='first'))
			constates[,i1] <- mind
		}
	}
	meanstates <- apply(constates, 1, mean)
	mcols(consensus.gr) <- meanstates

	# Variables
	nummod <- length(uni.hmm.list)
	file <- paste(file,"bed", sep=".")

	# Generate the colors
	colors <- state.colors[state.labels]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")

	# Write first line to file
	cat("browser hide all\n", file=file)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		uni.hmm <- uni.hmm.list[[imod]]
		uni.hmm.gr <- uni.hmm.grl[[imod]]
		priority <- 50 + imod
		cat(paste("track name=",uni.hmm$ID," description=\"univariate calls for ",names(uni.hmm.list)[imod],"\" visibility=1 itemRgb=On priority=",priority,"\n", sep=""), file=file, append=TRUE)

		# Change chromosome names from '1' to 'chr1' if necessary
		mask <- which(!grepl('chr', seqnames(uni.hmm.gr)))
		mcols(uni.hmm.gr)$chromosome <- as.character(seqnames(uni.hmm.gr))
		mcols(uni.hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(uni.hmm.gr)$chromosome[mask])
		mcols(uni.hmm.gr)$chromosome <- as.factor(mcols(uni.hmm.gr)$chromosome)

		# Collapse the calls
		collapsed.calls <- as.data.frame(uni.hmm.gr)[c('chromosome','start','end','state')]

		# Append to file
		itemRgb <- RGBs[as.character(collapsed.calls$state)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		write.table(format(df, scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

}


# ====================================
# Write signal tracks from binned data
# ====================================
binned2wiggle <- function(binned.data.list, file="view_me_in_genome_browser") {
	
	# Check user input
	if (is.null(names(binned.data.list))) {
		stop("Please name the list entries. The names will be used as track name for the Genome Browser file.")
	}

	# Variables
	file <- paste(file,"wig", sep=".")

	# Write to file
	cat("browser hide all\n", file=file)
	
	# Go through binned.data.list
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]
		names(binned.data) <- binned.data.names
		
		# Write track information
		name <- names(binned.data.list)[i1]
		binsize <- diff(binned.data$start[1:2])
		description <- paste("binned read counts ",binsize,"bp", sep="")
		cat(paste("track type=wiggle_0 name=\"",name,"\" description=\"",description,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",i1,"\n", sep=""), file=file, append=TRUE)
		# Write read data
		for (chrom in unique(binned.data$chrom)) {
			cat(paste("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n", sep=""), file=file, append=TRUE)
			write.table(binned.data$reads, file=file, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}

}


