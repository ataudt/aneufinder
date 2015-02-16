hmmList2GRangesList <- function(hmm.list, reduce=TRUE, numCPU=1, consensus=FALSE) {

	## Load models
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges
	cat('transforming to GRanges ...')
	if (numCPU > 1) {
		cl <- makeCluster(numCPU)
		registerDoParallel(cl)
		cfun <- function(...) { GRangesList(...) }
		hmm.grl <- foreach (hmm = hmm.list, .packages=c('aneufinder','GenomicRanges'), .combine='cfun', .multicombine=TRUE) %dopar% {
			hmm2GRanges(hmm, reduce=reduce)
		}
		if (consensus) {
			suppressMessages( consensus.gr <- disjoin(unlist(hmm.grl)) )
			constates <- foreach (hmm.gr = hmm.grl, .packages='GenomicRanges', .combine='cbind') %dopar% {
				splt <- split(hmm.gr, mcols(hmm.gr)$state)
				mind <- as.matrix(findOverlaps(consensus.gr, splt, select='first'))
			}
		}
		stopCluster(cl)
	} else {
		hmm.grl <- GRangesList()
		for (hmm in hmm.list) {
			hmm.grl[[length(hmm.grl)+1]] <- hmm2GRanges(hmm, reduce=reduce)
		}
		if (consensus) {
			consensus.gr <- disjoin(unlist(hmm.grl))
			constates <- matrix(NA, ncol=length(hmm.grl), nrow=length(consensus.gr))
			for (i1 in 1:length(hmm.grl)) {
				hmm.gr <- hmm.grl[[i1]]
				splt <- split(hmm.gr, mcols(hmm.gr)$state)
				mind <- as.matrix(findOverlaps(consensus.gr, splt, select='first'))
				constates[,i1] <- mind
			}
		}
	}
	cat(" done\n")
	if (consensus) {
		cat("calculating consensus template for all samples ...")
		meanstates <- apply(constates, 1, mean, na.rm=T)
		mcols(consensus.gr)$meanstate <- meanstates
		cat(" done\n")
		return(list(grl=hmm.grl, consensus=consensus.gr, constates=constates))
	} else {
		return(list(grl=hmm.grl))
	}


}

binned2GRanges <- function(binned.data, chrom.length.file=NULL, offset=0) {

	gr <- GenomicRanges::GRanges(
			seqnames = Rle(binned.data$chrom),
			ranges = IRanges(start=binned.data$start+offset, end=binned.data$end+offset),
			strand = Rle(strand("*"), nrow(binned.data)),
			reads = binned.data$reads
			)
	if (!is.null(chrom.length.file)) {
		# File with chromosome lengths (1-based)
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		seqlengths(gr) <- as.integer(chrom.lengths[names(seqlengths(gr))])
	}		
	return(gr)

}

hmm2GRanges <- function(hmm, reduce=TRUE) {

	### Check user input ###
	if (check.univariate.model(hmm)!=0 & check.multivariate.model(hmm)!=0) stop("argument 'hmm' expects a univariate or multivariate hmm object (type ?hmm for help)")
	if (check.logical(reduce)!=0) stop("argument 'reduce' expects TRUE or FALSE")

	### Create GRanges ###
	# Transfer coordinates
	gr <- GenomicRanges::GRanges(
			seqnames = Rle(hmm$coordinates$chrom),
			ranges = IRanges(start=hmm$coordinates$start, end=hmm$coordinates$end),
			strand = Rle(strand("*"), nrow(hmm$coordinates))
			)
	seqlengths(gr) <- hmm$seqlengths[names(seqlengths(gr))]
	# Reorder seqlevels
	gr <- keepSeqlevels(gr, names(hmm$seqlengths))

	if (reduce) {
		# Reduce state by state
		red.gr.list <- GenomicRanges::GRangesList()
		ustates <- unique(hmm$states)
		levels <- levels(hmm$states)
		for (state in ustates) {
			red.gr <- GenomicRanges::reduce(gr[hmm$states==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels),length(red.gr))
			if (class(hmm)==class.multivariate.hmm) {
				mcols(red.gr)$state.separate <- matrix(rep(strsplit(as.character(state), split=' ')[[1]], length(red.gr)), byrow=T, nrow=length(red.gr))
			}
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		# Merge and sort
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		remove(red.gr.list)
		return(red.gr)
	} else {
		mcols(gr)$reads <- hmm$reads
		mcols(gr)$posteriors <- hmm$posteriors
		mcols(gr)$state <- hmm$states
		if (class(hmm)==class.multivariate.hmm) {
			mcols(gr)$state.separate <- as.matrix(hmm$states.separate)
		}
		return(gr)
	}

}

bed2GRanges <- function(bedfile, chrom.length.file, skip=1, binsize=NULL) {

	# File with chromosome lengths (1-based)
	chrom.lengths.df <- read.table(chrom.length.file)
	chrom.lengths <- chrom.lengths.df[,2]
	names(chrom.lengths) <- chrom.lengths.df[,1]
	# File with reads, determine classes first for faster import (0-based)
	tab5rows <- read.table(bedfile, nrows=5, skip=skip)
	classes.in.bed <- sapply(tab5rows, class)
	classes <- rep("NULL",length(classes.in.bed))
	classes[1:4] <- classes.in.bed[1:4]
	data <- read.table(bedfile, colClasses=classes, skip=skip)
	# Convert to GRanges object
	data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]),
																	ranges=IRanges(start=data[,2]+1, end=data[,3]+1),	# +1 to match coordinate systems
																	strand=Rle(strand("*"), nrow(data)),
																	state=data[,4])
	seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	chroms.in.data <- seqlevels(data)

	## Inflate every range with bins
	if (!is.null(binsize)) {
		grl <- split(data, seqnames(data))
		inflated.data <- GRangesList()
		for (i1 in 1:length(grl)) {
			rgr <- ranges(grl[[i1]])
			widths <- width(rgr)
			numbins <- widths %/% binsize
			starts <- start(rgr)
			ends <- end(rgr)
			chroms <- seqnames(grl[[i1]])
			states <- mcols(grl[[i1]])$state

			# Create inflated vectors
			rle <- rle(1)
			rle$lengths <- numbins
			rle$values <- as.character(chroms)
			infchroms <- inverse.rle(rle)
			rle$values <- states
			infstates <- inverse.rle(rle)
			infstarts <- seq(starts[1], ends[length(ends)]-1, by=binsize)
			infends <- seq(starts[1]-1+binsize, ends[length(ends)]-2+binsize, by=binsize)

			inflated.chrom <- GenomicRanges::GRanges(seqnames=Rle(infchroms),
																								ranges=IRanges(start=infstarts, end=infends),
																								strand=Rle(strand('*'), sum(numbins)),
																								state=infstates)
			suppressWarnings( inflated.data[[i1]] <- inflated.chrom )
		}
		inflated.data <- unlist(inflated.data)
		seqlengths(inflated.data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		data <- inflated.data
	}
	mcols(data)$state <- factor(mcols(data)$state)

	return(data)

}
