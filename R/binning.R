bedGraph2binned <- function(bedGraphfile, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, gc.correction=TRUE, gc.correction.bsgenome, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bedGraphfile, format="bedGraph", chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, gc.correction=gc.correction, gc.correction.bsgenome=gc.correction.bsgenome, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bam2binned <- function(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, gc.correction=TRUE, gc.correction.bsgenome, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bamfile, format="bam", index=bamindex, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, gc.correction=gc.correction, gc.correction.bsgenome=gc.correction.bsgenome, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bed2binned <- function(bedfile, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, gc.correction=TRUE, gc.correction.bsgenome, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bedfile, format="bed", chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, gc.correction=gc.correction, gc.correction.bsgenome=gc.correction.bsgenome, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

align2binned <- function(file, format, index=file, chrom.length.file, outputfolder="binned_data", binsizes=NULL, reads.per.bin=10, numbins=NULL, chromosomes=NULL, gc.correction=TRUE, gc.correction.bsgenome, separate.chroms=FALSE, save.as.RData=TRUE) {

	## Uncomment this for use in debugging/developing
# 	format='bam'
# 	index=file
# 	binsizes=200000
# 	reads.per.bin=10
# 	numbins=NULL
# 	chromosomes=c(1:22,'X','Y')
# 	gc.correction=T
# 	separate.chroms=F
# 	save.as.RData=F
# 	library(BSgenome.Mmusculus.UCSC.mm10)
# 	gc.correction.bsgenome=BSgenome.Mmusculus.UCSC.mm10
# 	library(GenomicAlignments)

	## Check user input
	if (save.as.RData==FALSE) {
		separate.chroms=FALSE
	}
	if (gc.correction==TRUE) {
		check <- gc.correction.bsgenome	# trigger error if not defined
	}

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}

	### Read in the data
	## BED (0-based)
	if (format == "bed") {
		cat("Reading file",basename(file),"...")
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
		cat("Reading header of",basename(file),"...")
		file.header <- Rsamtools::scanBamHeader(file)[[1]]
		chrom.lengths <- file.header$targets
		chroms.in.data <- names(chrom.lengths)
	## BEDGraph (0-based)
	} else if (format == "bedGraph") {
		cat("Reading file",basename(file),"...")
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
	cat(" done\n")

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
 
	## Check if seqlengths of data and gc-correction are consistent
	if (gc.correction) {
		# Replace 1->chr1 if necessary
			chr.chroms2use <- chroms2use
			chr.chroms2use[!grepl('chr', chroms2use)] <- paste0('chr',chroms2use[!grepl('chr', chroms2use)])	
		# Compare
			compare <- chrom.lengths[chroms2use] == seqlengths(gc.correction.bsgenome)[chr.chroms2use]
			if (any(compare==FALSE)) {
				stop("Chromosome lengths differ between data and 'gc.correction.bsgenome'. Use the correct genome for option 'gc.correction.bsgenome'. You can also turn GC correction off by setting 'gc.correction=FALSE'.")
			}
	}

	## Determine binsize automatically
# 	if (!is.null(reads.per.bin)) {
		cat('Automatically determining binsizes...')
		gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use),
																	ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
		if (format=='bam') {
			autodata <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(what=c("pos"),which=range(gr),flag=scanBamFlag(isDuplicate=F)))
		} else if (format=='bed' | format=='bedGraph') {
			autodata <- data
		}
		numreadsperbp <- length(autodata) / sum(as.numeric(chrom.lengths[chroms2use]))
		## Pad binsizes and reads.per.bin with each others value
		binsizes.add <- round(reads.per.bin / numreadsperbp, -2)
		reads.per.bin.add <- round(binsizes * numreadsperbp, 2)
		binsizes <- c(binsizes, binsizes.add)
		reads.per.bin <- c(reads.per.bin.add, reads.per.bin)
		cat(' done\n')
# 	}

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
			cat("Binning into bin size",binsize,"with on average",readsperbin,"reads per bin\n")
		} else {
			numbin <- numbins[ibinsize]
			cat("Binning with number of bins",numbin,"\n")
		}

		### Iterate over all chromosomes
		binned.data <- GenomicRanges::GRanges()
		for (chromosome in chroms2use) {
			cat(chromosome,"                              \n")

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
				cat("initialize vectors...\r")
				chroms <- rep(chromosome,numbin)
				reads <- rep(0,numbin)
				start <- seq(from=1, by=binsize, length.out=numbin)
				end <- seq(from=binsize, by=binsize, length.out=numbin)
# 	 			end[length(end)] <- chrom.lengths[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired
			} else {
				## Initialize vectors
				cat("initialize vectors...\r")
				chroms <- rep(chromosome,numbin)
				reads <- rep(0,numbin)
				start <- round(seq(from=1, to=chrom.lengths[chromosome], length.out=numbin+1))
				end <- start[-1] - 1
				end[length(end)] <- end[length(end)] + 1
				start <- start[-length(start)]
			}

			## Create binned chromosome as GRanges object
			cat("creating GRanges container...            \r")
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbin),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbin)
							)
			seqlengths(i.binned.data) <- chrom.lengths[chromosome]

			if (format=="bam") {
				cat("reading reads from file...               \r")
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(what=c("pos"),which=range(i.binned.data),flag=scanBamFlag(isDuplicate=F)))
			}

			## Count overlaps
			if (format=="bam" | format=="bed") {
				cat("counting overlaps...                     \r")
				mreads <- GenomicRanges::countOverlaps(i.binned.data, data[seqnames(data)==chromosome & strand(data)=='-'])
				preads <- GenomicRanges::countOverlaps(i.binned.data, data[seqnames(data)==chromosome & strand(data)=='+'])
				reads <- mreads + preads
				
			} else if (format=="bedGraph") {
				cat("counting overlaps...                     \r")
				# Take the max value from all regions that fall into / overlap a given bin as read count
				midx <- as.matrix(findOverlaps(i.binned.data, data[seqnames(data)==chromosome]))
				reads <- rep(0,length(i.binned.data))
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
			
			### GC correction ###
			if (gc.correction) {
				cat("counting GC content...                   \r")

				library(Biostrings)
				# Correct seqnames 1->chr1 if necessary
				if (!grepl('chr',chromosome)) {
					chrom <- paste0('chr',chromosome)
				} else {
					chrom <- chromosome
				}

				## Calculating GC for whole bins
				view.chr <- Views(gc.correction.bsgenome[[chrom]], ranges(i.binned.data))
				freq <- alphabetFrequency(view.chr, as.prob = T, baseOnly=T)
				if (nrow(freq) > 1) {
					GC.bin <- rowSums(freq[, c("G","C")])
				} else {
					GC.bin <- sum(freq[, c("G","C")])
				}
				mcols(i.binned.data)$gc <- GC.bin
				
			}

			## Concatenate
			cat("concatenate...                           \r")
			if (is.null(numbins)) {
				mcols(i.binned.data)$reads <- reads
				mcols(i.binned.data)$mreads <- mreads
				mcols(i.binned.data)$preads <- preads
			} else {
				mcols(i.binned.data)$reads <- reads / (end - start)
				mcols(i.binned.data)$mreads <- mreads / (end - start)
				mcols(i.binned.data)$preads <- preads / (end - start)
			}

			if (separate.chroms==TRUE) {
				binned.data <- i.binned.data
				if (gc.correction) {
					cat("GC correction ...")
					binned.data$reads.gc <- binned.data$reads
					binned.data$mreads.gc <- binned.data$mreads
					binned.data$preads.gc <- binned.data$preads
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
						correction.factor <-  mean.reads.global / mean.reads.with.same.GC
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
					attr(binned.data, 'gc.correction') <- ggplt
					cat(" done\n")
				}
			
				if (save.as.RData==TRUE) {
					## Print to file
					if (is.null(numbins)) {
						filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F, trim=T),"_reads.per.bin_",readsperbin,"_",chromosome,"_.RData")
					} else {
						filename <- paste0(basename(file),"_numbin_",format(numbin, scientific=F, trim=T),"_",chromosome,"_.RData")
					}
					cat("save...                                  \r")
					save(binned.data, file=file.path(outputfolder,filename) )
				} else {
					cat("                                         \r")
					return(binned.data)
				}
			} else {
				binned.data <- suppressWarnings(BiocGenerics::append(binned.data, i.binned.data))
			}
			cat("                                         \r")

		}

		if (length(binned.data) == 0) {
			warning(paste0("The bin size of ",binsize," with reads per bin ",reads.per.bin," is larger than any of the chromosomes."))
			return(NULL)
		}

		if (gc.correction) {
			cat("GC correction ...")
			binned.data$reads.gc <- binned.data$reads
			binned.data$preads.gc <- binned.data$preads
			binned.data$mreads.gc <- binned.data$mreads
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
			attr(binned.data, 'gc.correction') <- ggplt
			cat(" done\n")
		}
			
		if (separate.chroms==FALSE) {
			if (save.as.RData==TRUE) {
				# Print to file
				if (is.null(numbins)) {
					filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F, trim=T),"_reads.per.bin_",readsperbin,"_.RData")
				} else {
					filename <- paste0(basename(file),"_numbin_",format(numbin, scientific=F, trim=T),"_.RData")
				}
				cat("Saving to file ...")
				save(binned.data, file=file.path(outputfolder,filename) )
				cat(" done\n")
			} else {
				return(binned.data)
			}
		}

	}

}
