# aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
# Copyright (C) 2015  Aaron Taudt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
#'
#' Convert aligned reads or signal values from various file formats into read counts in equidistant windows (bins). bam2binned() and bed2binned() use 'GenomicRanges::countOverlaps()' to calculate the read counts. Therefore, with small binsizes and large read lengths, one read fragment will often overlap more than one bin.
#' bedGraph2binned() sets the maximum signal value in a bin as value for that bin.
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
#' @param bamfile A file in BAM format.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @export
bam2binned <- function(bamfile, bamindex=bamfile, ID=basename(bamfile), pairedEndReads=FALSE, outputfolder="binned_data", binsizes=NULL, reads.per.bin=NULL, numbins=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, outputfolder.fragments=NULL, return.fragments=FALSE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bamfile, format="bam", ID=ID, bamindex=bamindex, pairedEndReads=pairedEndReads, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, save.as.RData=save.as.RData, calc.complexity=calc.complexity, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, call=call, outputfolder.fragments=outputfolder.fragments, return.fragments=return.fragments)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

#' @describeIn binning Bin reads in BED format
#' @inheritParams align2binned
#' @param bedfile A file in BED format.
#' @export
bed2binned <- function(bedfile, chrom.length.file, ID=basename(bedfile), outputfolder="binned_data", binsizes=NULL, reads.per.bin=NULL, numbins=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, outputfolder.fragments=NULL, return.fragments=FALSE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bedfile, format="bed", ID=ID, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, save.as.RData=save.as.RData, calc.complexity=calc.complexity, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, call=call, outputfolder.fragments=outputfolder.fragments, return.fragments=return.fragments)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

#' @describeIn binning Bin reads in bedGraph format
#' @inheritParams align2binned
#' @param bedGraphfile A file in bedGraph format.
#' @export
bedGraph2binned <- function(bedGraphfile, chrom.length.file, ID=basename(bedGraphfile), outputfolder="binned_data", binsizes=NULL, reads.per.bin=NULL, numbins=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, outputfolder.fragments=NULL, return.fragments=FALSE) {
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	binned.data <- align2binned(bedGraphfile, format="bedGraph", ID=ID, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, reads.per.bin=reads.per.bin, numbins=numbins, chromosomes=chromosomes, save.as.RData=save.as.RData, calc.complexity=calc.complexity, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, call=call, outputfolder.fragments=outputfolder.fragments, return.fragments=return.fragments)
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(binned.data)
}

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
#'
#' @param file A file with aligned reads.
#' @param format One of \code{c('bam', 'bed', 'bedGraph')}.
#' @param ID An identifier that will be used to identify the file throughout the workflow and in plotting.
#' @param bamindex Index file if \code{format='bam'} with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param chrom.length.file A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.
#' @param outputfolder Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.
#' @param binsizes A vector with integer values which will be used for the binning. If more than one value is given, output files will be produced for each bin size.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin size will be selected accordingly. Output files are produced for each value.
#' @param numbins Number of bins per chromosome. Each chromosome will have a different binsize! DO NOT USE UNLESS YOU KNOW WHAT YOU ARE DOING. Output files are produced for each value.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param save.as.RData If set to \code{FALSE}, no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case.
#' @param calc.complexity A logical indicating whether or not to estimate library complexity.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param remove.duplicate.reads A logical indicating whether or not duplicate reads should be removed.
#' @param call The \code{match.call()} of the parent function.
#' @param outputfolder.fragments Folder to which the read coordinates from the input file will be saved in \code{\link{GRanges}} format. These are the reads as used for the binning (duplicates removed, minimum mapping quality, etc.).
#' @param return.fragments If \code{TRUE} no binning is done and instead, read fragments from the input file are returned in \code{\link{GRanges}} format.
#' @return The function produces a \link{GRanges} object with one meta data column 'reads' that contains the read count. This binned data will be either written to file (\code{save.as.RData=FALSE}) or given as return value (\code{save.as.RData=FALSE}).
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairsFromBam readGAlignmentsFromBam first
#' @import preseqR
align2binned <- function(file, format, ID=basename(file), bamindex=file, pairedEndReads=FALSE, chrom.length.file, outputfolder="binned_data", binsizes=200000, reads.per.bin=NULL, numbins=NULL, chromosomes=NULL, save.as.RData=FALSE, calc.complexity=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, call=match.call(), outputfolder.fragments=NULL, return.fragments=FALSE) {

	## Check user input
	if (is.null(binsizes) & is.null(reads.per.bin)) {
		stop("Please specify either argument 'binsizes' or 'reads.per.bin'")
	}

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}
	if (!is.null(outputfolder.fragments)) {
		if (!file.exists(outputfolder.fragments)) {
			dir.create(outputfolder.fragments)
		}
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
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)))	# start+1 to go from [0,x) -> [1,x]
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		chroms.in.data <- seqlevels(data)
	## BAM (1-based)
	} else if (format == "bam") {
		## Check if bamindex exists
		bamindex.raw <- sub('\\.bai$', '', bamindex)
		bamindex <- paste0(bamindex.raw,'.bai')
		if (!file.exists(bamindex)) {
			bamindex.own <- Rsamtools::indexBam(file)
			warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
			bamindex <- bamindex.own
		}
		file.header <- Rsamtools::scanBamHeader(file)[[1]]
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
		gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
		if (calc.complexity || !remove.duplicate.reads) {
			if (pairedEndReads) {
				data <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq'))
				data <- GenomicAlignments::first(data)	# take only first mapping fragment of each pair
			} else {
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq'))
			}
		} else {
			if (pairedEndReads) {
				data <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
				data <- GenomicAlignments::first(data)	# take only first mapping fragment of each pair
			} else {
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
			}
		}
		## Filter by mapping quality
		if (!is.null(min.mapq)) {
			data <- data[mcols(data)$mapq >= min.mapq]
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
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)), signal=data[,4])	# start+1 to go from [0,x) -> [1,x]
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
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}
	diff <- setdiff(chromosomes, names(chrom.lengths))
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because no lengths could be found.'))
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	chroms2use <- intersect(chroms2use, names(chrom.lengths))
 
	if (calc.complexity) {
		message("  calculating complexity ...")
		downsample.sequence <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1)
		vm <- vector()
		k <- vector()
		multiplicity <- list()
		total.reads.sans.dup <- vector()
		total.reads.unique <- vector()
		total.reads <- vector()
		for (p in downsample.sequence) {
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
			dups['1'] <- length(down.data) - sum(dups*as.numeric(names(dups)))
			dups <- data.frame(multiplicity=as.numeric(names(dups)), frequency=dups)
			multiplicity[[as.character(p)]] <- dups

			total.reads.sans.dup[as.character(p)] <- dups[1,2]
			if (nrow(dups)>1) {
				total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)] + sum(dups[2:nrow(dups),2])
			} else {
				total.reads.unique[as.character(p)] <- total.reads.sans.dup[as.character(p)]
			}
			total.reads[as.character(p)] <- length(down.data)
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
	### Store the read fragments as GRanges ###
	fragments <- GRanges(seqnames=seqnames(data), ranges=IRanges(start=start(data), end=end(data)), strand=strand(data))
	if (!is.null(outputfolder.fragments)) {
		filename <- file.path(outputfolder.fragments,paste0(basename(file),"_fragments.RData"))
		save(fragments, file=filename)
	}
	if (return.fragments) {
		return(fragments)
	}

	### Loop over all binsizes ###
	## Pad binsizes and reads.per.bin with each others value
	numreadsperbp <- length(data) / sum(as.numeric(chrom.lengths[chroms2use]))
	binsizes.add <- round(reads.per.bin / numreadsperbp, -2)
	reads.per.bin.add <- round(binsizes * numreadsperbp, 2)
	binsizes <- c(binsizes, binsizes.add)
	reads.per.bin <- c(reads.per.bin.add, reads.per.bin)
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
		message("  binning genome ...", appendLF=F); ptm <- proc.time()
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

			suppressWarnings(
				binned.data[[chromosome]] <- i.binned.data
			)
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		binned.data <- unlist(binned.data)
		names(binned.data) <- NULL
		attr(binned.data, 'complexity.MM') <- NA
		attr(binned.data, 'complexity.MM.curve') <- NA
		attr(binned.data, 'complexity.preseqR') <- NA
		attr(binned.data, 'complexity.preseqR.curve') <- NA
		if (calc.complexity) {
			if (exists("complexity.MM.ggplt")) {
				attr(binned.data, 'complexity.MM.curve') <- complexity.MM.ggplt
			}
			if (exists("complexity.MM.fit")) {
				attr(binned.data, 'complexity.MM') <- complexity.MM
			}
			if (!is.null(complexity.preseqR.fit)) {
				attr(binned.data, 'complexity.preseqR') <- complexity.preseqR
				attr(binned.data, 'complexity.preseqR.curve') <- complexity.preseqR.ggplt
			}
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

		### Quality measures ###
		## Spikyness
		attr(binned.data, 'spikyness') <- qc.spikyness(binned.data$reads)
		## Shannon entropy
		attr(binned.data, 'shannon.entropy') <- qc.entropy(binned.data$reads)

		### ID ###
		attr(binned.data, 'ID') <- ID

		### Save or return the binned data ###
		if (save.as.RData==TRUE) {
			# Save to file
			if (is.null(numbins)) {
				filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F, trim=T),"_reads.per.bin_",readsperbin,"_.RData")
			} else {
				filename <- paste0(basename(file),"_numbin_",format(numbin, scientific=F, trim=T),"_.RData")
			}
			message("Saving to file ...", appendLF=F); ptm <- proc.time()
# 			attr(binned.data, 'call') <- call # do not store along with GRanges because it inflates disk usage
			save(binned.data, file=file.path(outputfolder,filename) )
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
# 			attr(binned.data, 'call') <- call
			return(binned.data)
		}

	}


}
