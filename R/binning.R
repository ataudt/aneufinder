bedGraph2binned <- function(bedGraphfile, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, GC.correction=TRUE, GC.correction.bsgenome, save.as.RData=TRUE, calc.complexity=TRUE, remove.duplicate.reads=TRUE, calc.spikyness=TRUE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bedGraphfile, format="bedGraph", pairedEndReads=FALSE, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, GC.correction=GC.correction, GC.correction.bsgenome=GC.correction.bsgenome, save.as.RData=save.as.RData, calc.complexity=calc.complexity, remove.duplicate.reads=remove.duplicate.reads, calc.spikyness=calc.spikyness, call=call)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

bam2binned <- function(bamfile, bamindex=bamfile, pairedEndReads=FALSE, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, GC.correction=TRUE, GC.correction.bsgenome, save.as.RData=TRUE, calc.complexity=TRUE, remove.duplicate.reads=TRUE, calc.spikyness=TRUE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bamfile, format="bam", index=bamindex, pairedEndReads=FALSE, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, GC.correction=GC.correction, GC.correction.bsgenome=GC.correction.bsgenome, save.as.RData=save.as.RData, calc.complexity=calc.complexity, remove.duplicate.reads=remove.duplicate.reads, calc.spikyness=calc.spikyness, call=call)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

bed2binned <- function(bedfile, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, GC.correction=TRUE, GC.correction.bsgenome, save.as.RData=TRUE, calc.complexity=TRUE, remove.duplicate.reads=TRUE, calc.spikyness=TRUE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bedfile, format="bed", pairedEndReads=FALSE, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, GC.correction=GC.correction, GC.correction.bsgenome=GC.correction.bsgenome, save.as.RData=save.as.RData, calc.complexity=calc.complexity, remove.duplicate.reads=remove.duplicate.reads, calc.spikyness=calc.spikyness, call=call)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

align2binned <- function(file, format, index=file, pairedEndReads=FALSE, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, GC.correction=TRUE, GC.correction.bsgenome, save.as.RData=TRUE, calc.complexity=TRUE, remove.duplicate.reads=TRUE, calc.spikyness=TRUE, call=match.call()) {

	## Uncomment this for use in debugging/developing
# 	format='bam'
# 	index=file
# 	outputfolder='binned_data'
# 	binsizes=200000
# 	reads.per.bin=10
# 	numbins=NULL
# 	chromosomes=c(1:22,'X','Y')
# 	GC.correction=T
# 	save.as.RData=F
# 	library(BSgenome.Mmusculus.UCSC.mm10)
# 	GC.correction.bsgenome=BSgenome.Mmusculus.UCSC.mm10
# 	library(BSgenome.Hsapiens.UCSC.hg19)
# 	GC.correction.bsgenome=BSgenome.Hsapiens.UCSC.hg19
# 	calc.complexity=T
# 	library(GenomicAlignments)
# 	library(ggplot2)
# 	library(GenomicRanges)

	## Check user input
	if (GC.correction==TRUE) {
		check <- GC.correction.bsgenome	# trigger error if not defined
	}

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}

	### Read in the data
	message("Reading file ",basename(file)," ...", appendLF=F); ptm <- proc.time()
	## BED (0-based)
	if (format == "bed") {
		# File with chromosome lengths (1-based)
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		# File with reads, determine classes first for faster import (0-based)
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[c(1:3,6)] <- classes.in.bed[c(1:3,6)]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]+1), strand=Rle(strand("*"), nrow(data)))	# +1 to match coordinate systems
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		chroms.in.data <- seqlevels(data)
	## BAM (1-based)
	} else if (format == "bam") {
		file.header <- Rsamtools::scanBamHeader(file)[[1]]
		chrom.lengths <- file.header$targets
		chroms.in.data <- names(chrom.lengths)
		if (is.null(chromosomes)) {
			chromosomes <- chroms.in.data
		}
		gr <- GenomicRanges::GRanges(seqnames=Rle(chromosomes),
																ranges=IRanges(start=rep(1, length(chromosomes)), end=chrom.lengths[chromosomes]))
		if (calc.complexity || !remove.duplicate.reads) {
			if (pairedEndReads) {
				data <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=index, param=ScanBamParam(which=range(gr)))
				data <- first(data)	# take only first mapping fragment of each pair
			} else {
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(which=range(gr)))
			}
		} else {
			if (pairedEndReads) {
				data <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=index, param=ScanBamParam(which=range(gr), flag=scanBamFlag(isDuplicate=F)))
				data <- first(data)	# take only first mapping fragment of each pair
			} else {
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(which=range(gr),flag=scanBamFlag(isDuplicate=F)))
			}
		}
	## BEDGraph (0-based)
	} else if (format == "bedGraph") {
		# File with chromosome lengths (1-based)
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		# File with reads, determine classes first for faster import
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:4] <- classes.in.bed[1:4]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]+1), strand=Rle(strand("*"), nrow(data)), signal=data[,4])	# +1 to match coordinate systems
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		chroms.in.data <- seqlevels(data)
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Select chromosomes to bin
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data'))
	}
	diff <- setdiff(chromosomes, names(chrom.lengths))
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because no lengths could be found'))
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	chroms2use <- intersect(chroms2use, names(chrom.lengths))
 
	## Check if seqlengths of data and GC.correction are consistent
	if (GC.correction) {
		# Replace 1->chr1 if necessary
			chr.chroms2use <- chroms2use
			chr.chroms2use[!grepl('chr', chroms2use)] <- paste0('chr',chroms2use[!grepl('chr', chroms2use)])	
		# Compare
			compare <- chrom.lengths[chroms2use] == seqlengths(GC.correction.bsgenome)[chr.chroms2use]
			if (any(compare==FALSE)) {
				stop("Chromosome lengths differ between data and 'GC.correction.bsgenome'. Use the correct genome for option 'GC.correction.bsgenome'. You can also turn GC correction off by setting 'GC.correction=FALSE'.")
			}
	}

	if (calc.complexity) {
		message("  calculating complexity ...")
		downsample.sequence <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1)
		vm <- vector()
		k <- vector()
		sum.unireads <- vector()
		sum.reads <- vector()
		for (i1 in 1:length(downsample.sequence)) {
			p <- downsample.sequence[i1]
			message("    p = ",p, appendLF=F)
			if (p != 1) {
				down.data <- data[sort(sample(1:length(data), size=p*length(data), replace=F))]
			} else {
				down.data <- data
			}
			sp <- start(down.data)[as.logical(strand(down.data)=='+')]
			sp1 <- c(sp[length(sp)], sp[-length(sp)])
			sm <- start(down.data)[as.logical(strand(down.data)=='-')]
			sm1 <- c(sm[length(sm)], sm[-length(sm)])
			sum.unireads[i1] <- length(which(sp!=sp1)) + length(which(sm!=sm1))
			sum.reads[i1] <- length(down.data)
		}
		message("")
		df <- data.frame(x=sum.reads, y=sum.unireads)
		vm.init <- quantile(df$y, 1)
		k.init <- quantile(df$x, .25)
		complexity.fit <- nls(y ~ vm * x/(k+x), data=df, start=list(vm=vm.init, k=k.init))
		x <- seq(from=0, to=5*max(sum.reads), length.out=100)
		df.fit <- data.frame(x=x, y=predict(complexity.fit, data.frame(x)))
		complexity.ggplt <- ggplot(df) + geom_point(aes(x=x, y=y), size=5) + geom_line(data=df.fit, mapping=aes(x=x, y=y)) + xlab('total number of reads') + ylab('uniquely mapping reads') + theme_bw()
		if (remove.duplicate.reads) {
			data <- c(data[strand(data)=='+'][sp!=sp1], data[strand(data)=='-'][sm!=sm1])
		}
	} else {
		if (remove.duplicate.reads) {
			sp <- start(data)[as.logical(strand(data)=='+')]
			sp1 <- c(sp[length(sp)], sp[-length(sp)])
			sm <- start(data)[as.logical(strand(data)=='-')]
			sm1 <- c(sm[length(sm)], sm[-length(sm)])
			data <- c(data[strand(data)=='+'][sp!=sp1], data[strand(data)=='-'][sm!=sm1])
		}
	}

	numreadsperbp <- length(data) / sum(as.numeric(chrom.lengths[chroms2use]))
	## Pad binsizes and reads.per.bin with each others value
	binsizes.add <- round(reads.per.bin / numreadsperbp, -2)
	reads.per.bin.add <- round(binsizes * numreadsperbp, 2)
	binsizes <- c(binsizes, binsizes.add)
	reads.per.bin <- c(reads.per.bin.add, reads.per.bin)

	### Do the loop for all binsizes
	if (is.null(numbins)) {
		length.binsizes <- length(binsizes)
	} else {
		length.binsizes <- length(numbins)
	}
	for (ibinsize in 1:length.binsizes) {
		if (is.null(numbins)) {
			binsize <- binsizes[ibinsize]
			readsperbin <- reads.per.bin[ibinsize]
			message("Binning into bin size ",binsize," with on average ",readsperbin," reads per bin")
		} else {
			numbin <- numbins[ibinsize]
			message("Binning with number of bins ",numbin)
		}

		### Iterate over all chromosomes
		if (GC.correction) {
			message("  binning genome with GC correction ...", appendLF=F); ptm <- proc.time()
		} else {
			message("  binning genome ...", appendLF=F); ptm <- proc.time()
		}
		binned.data <- GenomicRanges::GRangesList()
		for (chromosome in chroms2use) {

			if (is.null(numbins)) {
				## Check last incomplete bin
				incomplete.bin <- chrom.lengths[chromosome] %% binsize > 0
				if (incomplete.bin) {
					numbin <- floor(chrom.lengths[chromosome]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
				} else {
					numbin <- chrom.lengths[chromosome]/binsize
				}
				if (numbin == 0) {
					warning("chromosome ",chromosome," is smaller than binsize. Skipped.")
					next
				}
				## Initialize vectors
				chroms <- rep(chromosome,numbin)
				reads <- rep(0,numbin)
				start <- seq(from=1, by=binsize, length.out=numbin)
				end <- seq(from=binsize, by=binsize, length.out=numbin)
# 	 			end[length(end)] <- chrom.lengths[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired
			} else {
				## Initialize vectors
				chroms <- rep(chromosome,numbin)
				reads <- rep(0,numbin)
				start <- round(seq(from=1, to=chrom.lengths[chromosome], length.out=numbin+1))
				end <- start[-1] - 1
				end[length(end)] <- end[length(end)] + 1
				start <- start[-length(start)]
			}

			## Create binned chromosome as GRanges object
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbin),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbin)
							)
			seqlengths(i.binned.data) <- chrom.lengths[chromosome]

			### GC correction ###
			if (GC.correction) {
				# Correct seqnames 1->chr1 if necessary
				if (!grepl('chr',chromosome)) {
					chrom <- paste0('chr',chromosome)
				} else {
					chrom <- chromosome
				}
				## Calculating GC for whole bins
				view.chr <- Biostrings::Views(GC.correction.bsgenome[[chrom]], ranges(i.binned.data))
				freq <- Biostrings::alphabetFrequency(view.chr, as.prob = T, baseOnly=T)
				if (nrow(freq) > 1) {
					GC.bin <- rowSums(freq[, c("G","C")])
				} else {
					GC.bin <- sum(freq[, c("G","C")])
				}
				mcols(i.binned.data)$gc <- GC.bin
			}

			suppressWarnings(
				binned.data[[chromosome]] <- i.binned.data
			)
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		binned.data <- unlist(binned.data)
		names(binned.data) <- NULL
		if (calc.complexity) {
			attr(binned.data, 'complexity.ggplt') <- complexity.ggplt
			attr(binned.data, 'complexity.coefficients') <- coefficients(complexity.fit)
		}

		## Count overlaps
		message("  counting overlaps ...", appendLF=F); ptm <- proc.time()
		if (format=="bam" | format=="bed") {
			mreads <- GenomicRanges::countOverlaps(binned.data, data[strand(data)=='-'])
			preads <- GenomicRanges::countOverlaps(binned.data, data[strand(data)=='+'])
			reads <- mreads + preads
		} else if (format=="bedGraph") {
			# Take the max value from all regions that fall into / overlap a given bin as read count
			midx <- as.matrix(findOverlaps(binned.data, data))
			reads <- rep(0,length(binned.data))
			signal <- mcols(data)$signal
			rle <- rle(midx[,1])
			read.idx <- rle$values
			max.idx <- cumsum(rle$lengths)
			maxvalues <- rep(NA, length(read.idx))
			maxvalues[1] <- max(signal[midx[1:(max.idx[1]),2]])
			for (i1 in 2:length(read.idx)) {
				maxvalues[i1] <- max(signal[midx[(max.idx[i1-1]+1):(max.idx[i1]),2]])
			}
			reads[read.idx] <- maxvalues
			mreads <- rep(NA, length(reads))
			preads <- rep(NA, length(reads))
		}
		if (is.null(numbins)) {
			mcols(binned.data)$reads <- reads
			mcols(binned.data)$mreads <- mreads
			mcols(binned.data)$preads <- preads
		} else {
			mcols(binned.data)$reads <- reads / (end - start)
			mcols(binned.data)$mreads <- mreads / (end - start)
			mcols(binned.data)$preads <- preads / (end - start)
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		if (length(binned.data) == 0) {
			warning(paste0("The bin size of ",binsize," with reads per bin ",reads.per.bin," is larger than any of the chromosomes."))
			return(NULL)
		}

		if (GC.correction) {
			message("  GC correction ...", appendLF=F); ptm <- proc.time()
			binned.data$reads.gc <- binned.data$reads
			binned.data$preads.gc <- binned.data$preads
			binned.data$mreads.gc <- binned.data$mreads
			## Correction factors
			gc.categories <- seq(from=0, to=1, length=20)
			intervals.per.bin <- findInterval(binned.data$gc, gc.categories)
			intervals <- sort(unique(intervals.per.bin))
			mean.reads.global <- mean(binned.data$reads, trim=0.05)
			correction.factors <- NULL
			weights <- NULL
			for (interval in intervals) {
				mask <- intervals.per.bin==interval
				reads.with.same.GC <- binned.data$reads[mask]
				weights[as.character(gc.categories[interval])] <- length(reads.with.same.GC)
				mean.reads.with.same.GC <- mean(reads.with.same.GC, na.rm=T, trim=0.05)
				if (mean.reads.with.same.GC == 0) {
					correction.factor <- 0
				} else {
					correction.factor <-  mean.reads.global / mean.reads.with.same.GC
				}
				correction.factors[as.character(gc.categories[interval])] <- correction.factor
			}
			## Fit x^2 to correction.factors
			y <- correction.factors[-1][correction.factors[-1]<10]
			x <- as.numeric(names(y))
			w <- weights[-1][correction.factors[-1]<10]
			df <- data.frame(x,y,weight=w)
			fit <- lm(y ~ poly(x, 2, raw=T), data=df, weights=weight)
			fitted.correction.factors <- predict(fit, data.frame(x=gc.categories[intervals]))
			for (interval in intervals) {
				mask <- intervals.per.bin==interval
				correction.factor <- fitted.correction.factors[interval]
				binned.data$reads.gc[mask] <- binned.data$reads.gc[mask] * correction.factor
				binned.data$preads.gc[mask] <- binned.data$preads.gc[mask] * correction.factor
				binned.data$mreads.gc[mask] <- binned.data$mreads.gc[mask] * correction.factor
			}
			binned.data$reads.gc <- as.integer(round(binned.data$reads.gc))
			binned.data$preads.gc <- as.integer(round(binned.data$preads.gc))
			binned.data$mreads.gc <- as.integer(round(binned.data$mreads.gc))
			# Produce fit to check
			ggplt <- ggplot(df) + geom_point(aes(x=x, y=y, size=weight)) + geom_line(aes(x=x, y=y), data=data.frame(x=gc.categories[intervals], y=fitted.correction.factors)) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('correction factor')
			attr(binned.data, 'GC.correction.ggplt') <- ggplt
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}

		if (calc.spikyness) {
			reads <- binned.data$reads
			sum.reads <- sum(reads)
			spikyness <- vector()
			for (i1 in 1:10) {
				spikyness[i1] <- sum(abs(diff(sample(reads)))) / sum.reads
			}
			attr(binned.data, 'spikyness') <- mean(spikyness)
			if (GC.correction) {
				reads.gc <- binned.data$reads.gc
				sum.reads.gc <- sum(reads.gc)
				spikyness.gc <- vector()
				for (i1 in 1:10) {
					spikyness.gc[i1] <- sum(abs(diff(sample(reads.gc)))) / sum.reads.gc
				}
				attr(binned.data, 'spikyness.gc') <- mean(spikyness.gc)
			}
		}

		if (save.as.RData==TRUE) {
			# Print to file
			if (is.null(numbins)) {
				filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F, trim=T),"_reads.per.bin_",readsperbin,"_.RData")
			} else {
				filename <- paste0(basename(file),"_numbin_",format(numbin, scientific=F, trim=T),"_.RData")
			}
			message("Saving to file ...", appendLF=F); ptm <- proc.time()
			attr(binned.data, 'call') <- call
			save(binned.data, file=file.path(outputfolder,filename) )
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
			attr(binned.data, 'call') <- call
			return(binned.data)
		}

	}

}
