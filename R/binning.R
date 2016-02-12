

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows.
#'
#' Convert aligned reads or signal values from various file formats into read counts in equidistant windows (bins). bam2binned() and bed2binned() use 'GenomicRanges::countOverlaps()' to calculate the read counts. Therefore, with small binsizes and large read lengths, one read fragment will often overlap more than one bin.
#'
#' @name binning
#' @examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata/BB140820_I_002.bam", package="aneufinder")
#'## Bin the BAM file into bin size 200000bp
#'binned.data <- bam2binned(bamfile, binsize=200000, chromosomes=c(1:22,'X','Y'), save.as.RData=FALSE)
#' @author Aaron Taudt
NULL

#' @describeIn binning Bin reads in BAM format
#' @inheritParams align2binned
#' @inheritParams bam2GRanges
#' @export
bam2binned <- function(bamfile, bamindex=bamfile, ID=basename(bamfile), pairedEndReads=FALSE, max.fragment.width=1000, outputfolder.binned="binned_data", binsizes=NULL, reads.per.bin=NULL, stepsize=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, reads.store=FALSE, outputfolder.reads="data", reads.return=FALSE, reads.overwrite=FALSE, reads.only=FALSE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bamfile, format="bam", ID=ID, bamindex=bamindex, pairedEndReads=pairedEndReads, max.fragment.width=max.fragment.width, outputfolder.binned=outputfolder.binned, binsizes=binsizes, reads.per.bin=reads.per.bin, stepsize=stepsize, chromosomes=chromosomes, save.as.RData=save.as.RData, calc.complexity=calc.complexity, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, call=call, reads.store=reads.store, outputfolder.reads=outputfolder.reads, reads.return=reads.return, reads.overwrite=reads.overwrite, reads.only=reads.only)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

#' @describeIn binning Bin reads in BED format
#' @inheritParams align2binned
#' @param bedfile A file in BED format.
bed2binned <- function(bedfile, chrom.length.file, ID=basename(bedfile), outputfolder.binned="binned_data", binsizes=NULL, reads.per.bin=NULL, stepsize=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, reads.store=FALSE, outputfolder.reads="data", reads.return=FALSE, reads.overwrite=FALSE, reads.only=FALSE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bedfile, format="bed", ID=ID, chrom.length.file=chrom.length.file, outputfolder.binned=outputfolder.binned, binsizes=binsizes, reads.per.bin=reads.per.bin, stepsize=stepsize, chromosomes=chromosomes, save.as.RData=save.as.RData, calc.complexity=calc.complexity, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, call=call, reads.store=reads.store, outputfolder.reads=outputfolder.reads, reads.return=reads.return, reads.overwrite=reads.overwrite, reads.only=reads.only)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows.
#'
#' @param file A file with aligned reads.
#' @param format One of \code{c('bam', 'bed')}.
#' @param ID An identifier that will be used to identify the file throughout the workflow and in plotting.
#' @inheritParams bam2GRanges
#' @param chrom.length.file A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}).
#' @param outputfolder.binned Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.
#' @param binsizes An integer vector with bin sizes. If more than one value is given, output files will be produced for each bin size.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin size will be selected accordingly. Output files are produced for each value.
#' @param stepsize Fraction of the binsize that the sliding window is offset at each step. Example: If \code{stepsize=0.1} and \code{binsizes=c(200000,500000)}, the actual stepsize in basepairs is 20000 and 50000, respectively.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param save.as.RData If set to \code{FALSE}, no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case.
#' @param calc.complexity A logical indicating whether or not to estimate library complexity.
#' @param call The \code{match.call()} of the parent function.
#' @param reads.store If \code{TRUE} processed read fragments will be saved to file. Reads are processed according to \code{min.mapq} and \code{remove.duplicate.reads}. Paired end reads are coerced to single end fragments.
#' @param outputfolder.reads Folder to which the read fragments will be saved. If the specified folder does not exist, it will be created.
#' @param reads.return If \code{TRUE} no binning is done and instead, read fragments from the input file are returned in \code{\link{GRanges}} format.
#' @param reads.overwrite Whether or not an existing file with read fragments should be overwritten.
#' @param reads.only If \code{TRUE} only read fragments are stored and/or returned and no binning is done.
#' @return The function produces a \link{GRanges} object with one meta data column 'reads' that contains the read count. This binned data will be either written to file (\code{save.as.RData=FALSE}) or given as return value (\code{save.as.RData=FALSE}).
align2binned <- function(file, format, ID=basename(file), bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000, chrom.length.file, outputfolder.binned="binned_data", binsizes=200000, reads.per.bin=NULL, stepsize=NULL, save.as.RData=FALSE, calc.complexity=TRUE, call=match.call(), reads.store=FALSE, outputfolder.reads="data", reads.return=FALSE, reads.overwrite=FALSE, reads.only=FALSE) {

	## Check user input
	if (reads.return==FALSE & reads.only==FALSE) {
		if (is.null(binsizes) & is.null(reads.per.bin)) {
			stop("Please specify either argument 'binsizes' or 'reads.per.bin'")
		}
	}

	## Create outputfolder.binned if not exists
	if (!file.exists(outputfolder.binned) & save.as.RData==TRUE) {
		dir.create(outputfolder.binned)
	}

	### Read in the data
	message("Reading file ",basename(file)," ...", appendLF=FALSE); ptm <- proc.time()
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
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)))	# start+1 to go from [0,x) -> [1,x]
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	## BAM (1-based)
	} else if (format == "bam") {
		if (calc.complexity || !remove.duplicate.reads) {
			data <- suppressMessages( bam2GRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=FALSE, min.mapq=min.mapq, max.fragment.width=max.fragment.width) )
		} else {
			data <- suppressMessages( bam2GRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=TRUE, min.mapq=min.mapq, max.fragment.width=max.fragment.width) )
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Select chromosomes to bin
	if (is.null(chromosomes)) {
		chromosomes <- seqlevels(data)
	}
	chroms2use <- intersect(chromosomes, seqlevels(data))
 
	### Complexity estimation ###
	complexity <- NULL
	if (calc.complexity) {
		message("  calculating complexity ...", appendLF=FALSE); ptm <- proc.time()
		complexity <- suppressMessages( estimateComplexity(data) )
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	}
	if (remove.duplicate.reads) {
		sp <- start(data)[as.logical(strand(data)=='+')]
		sp1 <- c(sp[length(sp)], sp[-length(sp)])
		sm <- start(data)[as.logical(strand(data)=='-')]
		sm1 <- c(sm[length(sm)], sm[-length(sm)])
		data <- c(data[strand(data)=='+'][sp!=sp1], data[strand(data)=='-'][sm!=sm1])
	}

	### Return fragments if desired ###
	if (reads.store) {
		if (!file.exists(outputfolder.reads)) { dir.create(outputfolder.reads) }
		filename <- file.path(outputfolder.reads,paste0(ID,'.RData'))
		if (reads.overwrite | !file.exists(filename)) {
			save(data, file=filename)
		}
	}
	if (reads.return) {
		return(data)
	}
	if (reads.only) {
		return()
	}

	### Coverage and percentage of genome covered ###
	genome.length <- sum(as.numeric(seqlengths(data)))
	data.strand <- data
	strand(data.strand) <- '*'
	coverage <- sum(as.numeric(width(data.strand))) / genome.length
	genome.covered <- sum(as.numeric(width(reduce(data.strand)))) / genome.length
	## Per chromosome
	coverage.per.chrom <- vector()
	genome.covered.per.chrom <- vector()
	for (chr in chroms2use) {
		data.strand.chr <- data.strand[seqnames(data.strand)==chr]
		coverage.per.chrom[chr] <- sum(as.numeric(width(data.strand.chr))) / seqlengths(data.strand)[chr]
		genome.covered.per.chrom[chr] <- sum(as.numeric(width(reduce(data.strand.chr)))) / seqlengths(data.strand)[chr]
	}
	coverage <- list(coverage=coverage, genome.covered=genome.covered, coverage.per.chrom=coverage.per.chrom, genome.covered.per.chrom=genome.covered.per.chrom)


	## Pad binsizes and reads.per.bin with each others value
	numcountsperbp <- length(data) / sum(as.numeric(seqlengths(data)))
	binsizes.add <- round(reads.per.bin / numcountsperbp, -2)
	reads.per.bin.add <- round(binsizes * numcountsperbp, 2)
	binsizes <- c(binsizes, binsizes.add)
	reads.per.bin <- c(reads.per.bin.add, reads.per.bin)
	length.binsizes <- length(binsizes)

	### Loop over all binsizes ###
	data.plus <- data[strand(data)=='+']
	data.minus <- data[strand(data)=='-']
	skipped.chroms <- character()
	for (ibinsize in 1:length.binsizes) {
		binsize <- binsizes[ibinsize]
		readsperbin <- reads.per.bin[ibinsize]
		message("Binning into bin size ",binsize," with on average ",readsperbin," reads per bin")

		### Loop over chromosomes ###
		message("  binning genome ...", appendLF=FALSE); ptm <- proc.time()
		binned.data <- GenomicRanges::GRangesList()
		for (chromosome in chroms2use) {

			## Check last incomplete bin
			incomplete.bin <- seqlengths(data)[chromosome] %% binsize > 0
			if (incomplete.bin) {
				numbin <- floor(seqlengths(data)[chromosome]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
			} else {
				numbin <- seqlengths(data)[chromosome]/binsize
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
# 			end[length(end)] <- seqlengths(data)[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired

			## Create binned chromosome as GRanges object
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbin),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbin)
							)
			seqlengths(i.binned.data) <- seqlengths(data)[chromosome]

			suppressWarnings(
				binned.data[[chromosome]] <- i.binned.data
			)

		} ### end loop chromosomes ###
		if (length(skipped.chroms)>0) {
			warning("Chromosomes ", paste0(skipped.chroms, collapse=', '), " are smaller than binsize. Skipped.")
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		binned.data <- unlist(binned.data)
		names(binned.data) <- NULL

		### Loop over offsets ###
		countmatrices <- list()
		if (is.null(stepsize)) {
			offsets <- 0
		} else {
			offsets <- seq(from=0, to=binsize, by=as.integer(stepsize*binsize))
		}
		for (ioff in offsets) {
			## Count overlaps
			message(paste0("  counting overlaps for offset ",ioff," ..."), appendLF=FALSE); ptm <- proc.time()
			binned.data.shift <- suppressWarnings( shift(binned.data, shift=ioff) )
			if (format=="bam" | format=="bed") {
				mcounts <- suppressWarnings( GenomicRanges::countOverlaps(binned.data.shift, data.minus) )
				pcounts <- suppressWarnings( GenomicRanges::countOverlaps(binned.data.shift, data.plus) )
				counts <- mcounts + pcounts
			}
			countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
			colnames(countmatrix) <- c('counts','mcounts','pcounts')
			countmatrices[[as.character(ioff)]] <- countmatrix
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		mcols(binned.data) <- as(countmatrices[['0']],'DataFrame')
		attr(binned.data,'offset.counts') <- countmatrices

		if (length(binned.data) == 0) {
			warning(paste0("The bin size of ",binsize," with reads per bin ",reads.per.bin," is larger than any of the chromosomes."))
			return(NULL)
		}

		### Quality measures ###
		qualityInfo <- list(complexity=complexity, coverage=coverage, spikyness=qc.spikyness(binned.data$counts), shannon.entropy=qc.entropy(binned.data$counts))
		attr(binned.data, 'qualityInfo') <- qualityInfo
		attr(binned.data, 'min.mapq') <- min.mapq

		### ID ###
		attr(binned.data, 'ID') <- ID

		### Save or return the binned data ###
		if (save.as.RData==TRUE) {
			# Save to file
			filename <- paste0(ID,"_binsize_",format(binsize, scientific=TRUE, trim=TRUE),"_reads.per.bin_",readsperbin,"_.RData")
			message("Saving to file ...", appendLF=FALSE); ptm <- proc.time()
# 			attr(binned.data, 'call') <- call # do not store along with GRanges because it inflates disk usage
			save(binned.data, file=file.path(outputfolder.binned,filename) )
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
# 			attr(binned.data, 'call') <- call
			return(binned.data)
		}

	} ### end loop binsizes ###


}


#' Estimate library complexity
#'
#' Estimate library complexity using a very simple "Michaelis-Menten" approach and the sophisticated approach from the \pkg{\link{preseqR}} package.
#'
#' @param reads A \code{\link{GRanges}} object with read fragments. NOTE: Complexity estimation relies on duplicate reads and therefore the duplicates have to be present in the input.
#' @import preseqR
estimateComplexity <- function(reads) {
	message("Calculating complexity")
	downsample.sequence <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1) # downsampling sequence for MM approach
	vm <- vector()
	k <- vector()
	multiplicity <- list()
	total.reads.sans.dup <- vector()
	total.reads.unique <- vector()
	total.reads <- vector()
	for (p in downsample.sequence) {
		message("    p = ",p, appendLF=FALSE)
		if (p != 1) {
			down.reads <- reads[sort(sample(1:length(reads), size=p*length(reads), replace=FALSE))]
		} else {
			down.reads <- reads
		}
		sp <- start(down.reads)[as.logical(strand(down.reads)=='+')]
		sp1 <- c(sp[length(sp)], sp[-length(sp)])
		sm <- start(down.reads)[as.logical(strand(down.reads)=='-')]
		sm1 <- c(sm[length(sm)], sm[-length(sm)])

		if (length(sp)==1) {
			rlep <- rle(FALSE)
		} else {
			rlep <- rle(sp==sp1)
		}
		if (length(sm)==1) {
			rlem <- rle(FALSE)
		} else {
			rlem <- rle(sm==sm1)
		}
		tab.p <- table(rlep$lengths[rlep$values])	# table of number of duplicates
		names(tab.p) <- as.numeric(names(tab.p)) + 1
		tab.m <- table(rlem$lengths[rlem$values])
		names(tab.m) <- as.numeric(names(tab.m)) + 1
		multiplicities <- sort(as.numeric(c(1, union(names(tab.p), names(tab.m)))))
		m <- matrix(0, nrow=length(multiplicities), ncol=2)
		rownames(m) <- multiplicities
		m[names(tab.p),1] <- tab.p
		m[names(tab.m),2] <- tab.m
		dups <- apply(m, 1, sum)
		dups['1'] <- length(down.reads) - sum(dups*as.numeric(names(dups)))
		dups <- data.frame(multiplicity=as.numeric(names(dups)), frequency=dups)
		multiplicity[[as.character(p)]] <- dups

		total.reads.sans.dup[as.character(p)] <- dups[1,2]
		if (nrow(dups)>1) {
			total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)] + sum(dups[2:nrow(dups),2])
		} else {
			total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)]
		}
		total.reads[as.character(p)] <- length(down.reads)
	}
	message("")
	df <- data.frame(x=total.reads, y=total.reads.unique)

	## Complexity estimation with preseqR
	complexity.preseqR.fit <- preseqR::preseqR.rfa.curve(multiplicity[['1']])
	if (is.null(complexity.preseqR.fit)) {
		complexity.preseqR <- NA
		complexity.preseqR.ggplt <- NA
		warning("Complexity estimation with preseqR failed.")
	} else {
		complexity.preseqR <- as.numeric(preseqR::preseqR.rfa.estimate(complexity.preseqR.fit$continued.fraction, 1e100) + total.reads.unique['1'])
		complexity.preseqR.ggplt <- ggplot(as.data.frame(complexity.preseqR.fit$estimates)) + geom_line(aes_string(x='sample.size', y='yield.estimate')) + geom_point(data=df, aes_string(x='x', y='y'), size=5) + ggtitle('preseqR complexity estimation') + theme_bw() + xlab('total number of reads') + ylab('unique reads')
	}

	## Complexity estimation with Michaelis-Menten
	complexity.MM <- NA
	complexity.MM.ggplt <- NA
	vm.init <- quantile(df$y, 1)
	k.init <- quantile(df$x, .25)
	tC <- tryCatch({
		complexity.MM.fit <- nls(y ~ vm * x/(k+x), data=df, start=list(vm=vm.init, k=k.init))
		complexity.MM <- as.numeric(coefficients(complexity.MM.fit)[1])
		max.x <- max(0.9 * coefficients(complexity.MM.fit)['k'] / (1-0.9), max(total.reads))
		x <- seq(from=0, to=max.x, length.out=1000)
		df.fit <- data.frame(x=x, y=predict(complexity.MM.fit, data.frame(x)))
		complexity.MM.ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y'), size=5) + geom_line(data=df.fit, mapping=aes_string(x='x', y='y')) + xlab('total number of reads') + ylab('unique reads') + theme_bw() + ggtitle('Michaelis-Menten complexity estimation')
	}, error = function(err) {
		warning("Complexity estimation with Michaelis-Menten failed.")
	})

	rl <- list(preseqR=complexity.preseqR, preseqR.ggplt=complexity.preseqR.ggplt, MM=complexity.MM, MM.ggplt=complexity.MM.ggplt)
	return(rl)
}


#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param bamfile A file in BAM format.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairsFromBam readGAlignmentsFromBam first
bam2GRanges <- function(bamfile, bamindex=bamfile, chromosomes=NULL, pairedEndReads=FALSE, remove.duplicate.reads=FALSE, min.mapq=10, max.fragment.width=1000) {

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
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=FALSE)))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=FALSE)))
		}
	}
	if (length(data.raw) == 0) {
		if (pairedEndReads) {
			stop(paste0("No reads imported. Does your file really contain paired end reads? Try with 'pairedEndReads=FALSE'"))
		}
		stop(paste0('No reads imported! Check your BAM-file ', bamfile))
	}
	## Filter by mapping quality
	if (pairedEndReads) {
		message("Converting to GRanges ...", appendLF=FALSE); ptm <- proc.time()
		data <- as(data.raw, 'GRanges') # treat as one fragment
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		message("Filtering reads ...", appendLF=FALSE); ptm <- proc.time()
		if (!is.null(min.mapq)) {
			mapq.first <- mcols(GenomicAlignments::first(data.raw))$mapq
			mapq.last <- mcols(GenomicAlignments::last(data.raw))$mapq
			mapq.mask <- mapq.first >= min.mapq & mapq.last >= min.mapq
			if (any(is.na(mapq.mask))) {
				warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
			}
			data <- data[which(mapq.mask)]
			# Filter out too long fragments
			data <- data[width(data)<=max.fragment.width]
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	} else {
		message("Converting to GRanges ...", appendLF=FALSE); ptm <- proc.time()
		data <- as(data.raw, 'GRanges')
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		message("Filtering reads ...", appendLF=FALSE); ptm <- proc.time()
		if (!is.null(min.mapq)) {
			if (any(is.na(mcols(data)$mapq))) {
				warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
				mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
			}
			data <- data[mcols(data)$mapq >= min.mapq]
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	}
	return(data)

}

