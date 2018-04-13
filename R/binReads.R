

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed(.gz) format into read counts in equidistant windows.
#'
#' Convert aligned reads from .bam or .bed(.gz) files into read counts in equidistant windows (bins). This function uses \code{GenomicRanges::countOverlaps} to calculate the read counts.
#'
#' @param file A file with aligned reads. Alternatively a \code{\link{GRanges-class}} with aligned reads.
#' @param ID An identifier that will be used to identify the file throughout the workflow and in plotting.
#' @inheritParams bam2GRanges
#' @inheritParams bed2GRanges
#' @param outputfolder.binned Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.
#' @param binsizes An integer vector with bin sizes. If more than one value is given, output files will be produced for each bin size.
#' @param stepsizes A vector of step sizes the same length as \code{binsizes}. Only used for \code{method="HMM"}.
#' @param bins A named \code{list} with \code{\link{GRanges-class}} containing precalculated bins produced by \code{\link{fixedWidthBins}} or \code{\link{variableWidthBins}}. Names must correspond to the binsize.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin size will be selected accordingly. Output files are produced for each value.
#' @param reads.per.step Approximate number of desired reads per step.
#' @param variable.width.reference A BAM file that is used as reference to produce variable width bins. See \code{\link{variableWidthBins}} for details.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param save.as.RData If set to \code{FALSE}, no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case.
#' @param calc.complexity A logical indicating whether or not to estimate library complexity.
#' @param call The \code{match.call()} of the parent function.
#' @param reads.store If \code{TRUE} processed read fragments will be saved to file. Reads are processed according to \code{min.mapq} and \code{remove.duplicate.reads}. Paired end reads are coerced to single end fragments. Will be ignored if \code{use.bamsignals=TRUE}.
#' @param outputfolder.reads Folder to which the read fragments will be saved. If the specified folder does not exist, it will be created.
#' @param reads.return If \code{TRUE} no binning is done and instead, read fragments from the input file are returned in \code{\link{GRanges-class}} format.
#' @param reads.overwrite Whether or not an existing file with read fragments should be overwritten.
#' @param reads.only If \code{TRUE} only read fragments are stored and/or returned and no binning is done.
#' @param use.bamsignals If \code{TRUE} the \pkg{\link[bamsignals]{bamsignals}} package will be used for binning. This gives a tremendous performance increase for the binning step. \code{reads.store} and \code{calc.complexity} will be set to \code{FALSE} in this case.
#' @return The function produces a \code{list()} of \code{\link{GRanges-class}} or \code{\link{GRangesList}} objects with meta data columns 'counts', 'mcounts', 'pcounts' that contain the total, minus and plus read count. This binned data will be either written to file (\code{save.as.RData=FALSE}) or given as return value (\code{save.as.RData=FALSE}).
#' @seealso binning
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom bamsignals bamCount
#' @export
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the BED file into bin size 1Mb
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'print(binned)
#'
binReads <- function(file, assembly, ID=basename(file), bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000, blacklist=NULL, outputfolder.binned="binned_data", binsizes=1e6, stepsizes=NULL, reads.per.bin=NULL, reads.per.step=NULL, bins=NULL, variable.width.reference=NULL, save.as.RData=FALSE, calc.complexity=TRUE, call=match.call(), reads.store=FALSE, outputfolder.reads="data", reads.return=FALSE, reads.overwrite=FALSE, reads.only=FALSE, use.bamsignals=FALSE) {

	## Determine format
	if (is.character(file)) {
			file.clean <- sub('\\.gz$','', file)
			format <- rev(strsplit(file.clean, '\\.')[[1]])[1]
	} else if (class(file)=='GRanges') {
			format <- 'GRanges'
	}

	if (format=='bed') {
			temp <- assembly # trigger error if not defined
	}

	## Variables for bamsignals
	paired.end <- 'ignore'
	if (pairedEndReads) {
			paired.end <- 'filter'
	}
	filteredFlag <- -1
	if (remove.duplicate.reads) {
	    filteredFlag <- 1024
	}
	if (use.bamsignals) {
	    reads.store <- FALSE
	    calc.complexity <- FALSE
	}

	## Check user input
	if (reads.return==FALSE & reads.only==FALSE) {
		if (is.null(binsizes) & is.null(reads.per.bin) & is.null(bins)) {
			stop("Please specify either argument 'binsizes' or 'reads.per.bin'")
		}
	}
	if (!is.null(stepsizes)) {
	    if (any(binsizes < stepsizes)) {
	        stop("'stepsizes' must be smaller/equal than 'binsizes'")
	    }
	}

	## Create outputfolder.binned if not exists
	if (!file.exists(outputfolder.binned) & save.as.RData==TRUE) {
		dir.create(outputfolder.binned)
	}

	### Read in the data
	data <- NULL
	if (format == "bed") {
		## BED (0-based)
		if (calc.complexity || !remove.duplicate.reads) {
			data <- bed2GRanges(file, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=FALSE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
		} else {
			data <- bed2GRanges(file, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=TRUE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
		}
		chrom.lengths <- seqlengths(data)
	} else if (format == "bam") {
		## BAM (1-based)
		if (use.bamsignals) {
			## Check if bamindex exists
			bamindex.raw <- sub('\\.bai$', '', bamindex)
			bamindex <- paste0(bamindex.raw,'.bai')
			if (!file.exists(bamindex)) {
					ptm <- startTimedMessage("Making bam-index file ...")
					bamindex.own <- Rsamtools::indexBam(file)
					# warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
					bamindex <- bamindex.own
					stopTimedMessage(ptm)
			}
			ptm <- startTimedMessage(paste0("Reading header from ", file, " ..."))
			chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(file))
			stopTimedMessage(ptm)
		} else {
			if (calc.complexity || !remove.duplicate.reads) {
				data <- bam2GRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=FALSE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
			} else {
				data <- bam2GRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=TRUE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
			}
			chrom.lengths <- seqlengths(data)
		}
	} else if (format == "GRanges") {
		## GRanges (1-based)
		data <- file
		err <- tryCatch({
			!is.character(ID)
		}, error = function(err) {
			TRUE
		})
		if (err) {
			ID <- 'GRanges'
		}
		chrom.lengths <- seqlengths(data)
	}

	## Select chromosomes to bin
	if (is.null(chromosomes)) {
			chromosomes <- names(chrom.lengths)
	}
	chroms2use <- intersect(chromosomes, names(chrom.lengths))
	skipped.chroms <- setdiff(chromosomes, chroms2use)
	## Stop if none of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring, '. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	if (length(skipped.chroms)>0) {
			warning("Could not find chromosomes ", paste0(skipped.chroms, collapse=', '), ".")
	}

	if (!is.null(data)) {
		### Select only desired chromosomes ###
		ptm <- startTimedMessage("Subsetting specified chromosomes ...")
		data <- data[seqnames(data) %in% chroms2use]
		data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
		## Drop seqlevels where seqlength is NA
		na.seqlevels <- seqlevels(data)[is.na(seqlengths(data))]
		data <- data[seqnames(data) %in% seqlevels(data)[!is.na(seqlengths(data))]]
		data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
		if (length(na.seqlevels) > 0) {
				warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
		}
		stopTimedMessage(ptm)
  
		### Complexity estimation ###
		complexity <- c(MM=NA)
		if (calc.complexity) {
			ptm <- startTimedMessage("Calculating complexity ...")
			complexity <- suppressMessages( estimateComplexity(data)[[1]] )
			stopTimedMessage(ptm)
		}
		if (remove.duplicate.reads) {
			ptm <- startTimedMessage("Removing duplicate reads ...")
			sp <- start(data)[as.logical(strand(data)=='+')]
			sp1 <- c(sp[length(sp)], sp[-length(sp)])
			sm <- start(data)[as.logical(strand(data)=='-')]
			sm1 <- c(sm[length(sm)], sm[-length(sm)])
			data <- c(data[strand(data)=='+'][sp!=sp1], data[strand(data)=='-'][sm!=sm1])
			stopTimedMessage(ptm)
		}

		### Return fragments if desired ###
		if (reads.store) {
			if (!file.exists(outputfolder.reads)) { dir.create(outputfolder.reads) }
			filename <- file.path(outputfolder.reads,paste0(ID,'.RData'))
			if (reads.overwrite | !file.exists(filename)) {
				ptm <- startTimedMessage(paste0("Saving fragments to ", filename, " ..."))
				save(data, file=filename)
				stopTimedMessage(ptm)
			}
		}
		if (reads.return) {
			return(data)
		}
		if (reads.only) {
		  return()
		}

		### Coverage and percentage of genome covered ###
		ptm <- startTimedMessage("Calculating coverage ...")
		genome.length <- sum(as.numeric(seqlengths(data)))
		data.strand <- data
		strand(data.strand) <- '*'
		coverage <- sum(as.numeric(width(data.strand))) / genome.length
		genome.covered <- sum(as.numeric(width(reduce(data.strand)))) / genome.length
		## Per chromosome
		coverage.per.chrom <- numeric()
		genome.covered.per.chrom <- numeric()
		for (chr in chroms2use) {
			data.strand.chr <- data.strand[seqnames(data.strand)==chr]
			coverage.per.chrom[chr] <- sum(as.numeric(width(data.strand.chr))) / seqlengths(data.strand)[chr]
			genome.covered.per.chrom[chr] <- sum(as.numeric(width(reduce(data.strand.chr)))) / seqlengths(data.strand)[chr]
		}
		coverage <- list(coverage=coverage, genome.covered=genome.covered, coverage.per.chrom=coverage.per.chrom, genome.covered.per.chrom=genome.covered.per.chrom)
		stopTimedMessage(ptm)

	} else {
		complexity <- c(MM=NA)
		coverage <- NA
	}

	### Make variable width bins ###
	if (!is.null(variable.width.reference)) {
			message("Making variable width bins:")
			if (is.character(variable.width.reference)) {
					variable.width.reference.clean <- sub('\\.gz$','', variable.width.reference)
					vformat <- rev(strsplit(variable.width.reference.clean, '\\.')[[1]])[1]
			} else if (class(variable.width.reference)=='GRanges') {
					vformat <- 'GRanges'
			}
			if (vformat == 'bam') {
					refreads <- bam2GRanges(variable.width.reference, bamindex=variable.width.reference, chromosomes=chroms2use, pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
			} else if (vformat == 'bed') {
					refreads <- bed2GRanges(variable.width.reference, assembly=assembly, chromosomes=chroms2use, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
			}
			bins.binsize <- variableWidthBins(refreads, binsizes=binsizes, stepsizes=stepsizes, chromosomes=chroms2use)
			message("Finished making variable width bins.")
	}
			
	### Make fixed width bins ###
	if (is.null(variable.width.reference)) {
			bins.binsize <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes, stepsizes=stepsizes)
	}

	### Make reads.per.bin bins ###
	bins.rpb <- NULL
	if (!is.null(data)) {
			numcountsperbp <- length(data) / sum(as.numeric(seqlengths(data)))
			binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
			stepsizes.rpb <- round(reads.per.step / numcountsperbp, -2)
			bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb, stepsizes=stepsizes.rpb)
	} else if (use.bamsignals & !is.null(reads.per.bin)) {
			ptm <- startTimedMessage("Parsing bamfile to determine binsize for reads.per.bin option ...")
			bins.helper <- suppressMessages( fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=1e6)[[1]] )
			counts <- bamsignals::bamCount(file, bins.helper, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE)
			stopTimedMessage(ptm)
			numcountsperbp <- sum(as.numeric(counts)) / sum(as.numeric(chrom.lengths[chroms2use]))
			binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
			stepsizes.rpb <- round(reads.per.step / numcountsperbp, -2)
			bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb, stepsizes=stepsizes.rpb)
	}
	
	### Combine in bins.list ###
	bins.list <- c(bins, bins.binsize, bins.rpb)

	### Loop over all binsizes ###
	if (!use.bamsignals | format == 'bed') {
			ptm <- startTimedMessage("Splitting into strands ...")
			data.plus <- data[strand(data)=='+']
			data.minus <- data[strand(data)=='-']
			data.star <- data[strand(data)=='*']
			ptm <- stopTimedMessage(ptm)
	}
	for (ibinsize in 1:length(bins.list)) {
			binsize <- as.numeric(sub('binsize_', '', sub('_stepsize.*', '', names(bins.list)[ibinsize])))
			ptm <- startTimedMessage("Counting overlaps for ", names(bins.list)[ibinsize], " ...")
			bins.steplist <- bins.list[[ibinsize]]
			if (class(bins.list[[ibinsize]]) == 'GRanges') {
			  bins.steplist <- GRangesList('0'=bins.steplist)
			}
			bins.steplist.counts <- GRangesList()
			for (istep in 1:length(bins.steplist)) {
			    bins <- bins.steplist[[istep]]
    			if (length(bins) == 0) {
    				stop(paste0("The bin size of ",binsize," with reads per bin ",reads.per.bin," is larger than any of the chromosomes."))
    			}

    			if (format == 'bam' & use.bamsignals) {
    					counts.ss <- bamsignals::bamCount(file, bins, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE, ss=TRUE, filteredFlag=filteredFlag)
    					pcounts <- counts.ss[1,]
    					mcounts <- counts.ss[2,]
    					counts <- pcounts + mcounts
    					readsperbin <- round(sum(as.numeric(counts)) / length(counts), 2)
    					
    			} else {
    					readsperbin <- round(length(data) / sum(as.numeric(seqlengths(data))) * binsize, 2)
    					scounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.star) )
    					mcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.minus) )
    					pcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.plus) )
    					counts <- mcounts + pcounts + scounts
    	
    			}
    			countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
    			colnames(countmatrix) <- c('counts','mcounts','pcounts')
    			mcols(bins) <- as(countmatrix, 'DataFrame')
    			bins.steplist.counts[[istep]] <- bins
			}
					
			if (class(bins.list[[ibinsize]]) == 'GRanges') {
  			bins <- bins.steplist.counts[[1]]
			} else {
  			bins <- bins.steplist.counts
			}
			
			### Quality measures ###
			qualityInfo <- list(complexity=complexity, coverage=coverage, spikiness=qc.spikiness(bins$counts), entropy=qc.entropy(bins$counts))
			attr(bins, 'qualityInfo') <- qualityInfo
			attr(bins, 'min.mapq') <- min.mapq

			### ID ###
			attr(bins, 'ID') <- ID
			stopTimedMessage(ptm)

			### Save or return the bins ###
			if (save.as.RData==TRUE) {
				# Save to file
				filename <- paste0(ID,"_", names(bins.list)[ibinsize],"_reads.per.bin_",readsperbin,"_.RData")
				ptm <- startTimedMessage("Saving to file ...")
				save(bins, file=file.path(outputfolder.binned,filename) )
				stopTimedMessage(ptm)
			} else {
				bins.list[[ibinsize]] <- bins
			}

	} ### end loop binsizes ###

	if (!save.as.RData) {
		return(bins.list)
	}


}


#' Estimate library complexity
#'
#' Estimate library complexity using a very simple "Michaelis-Menten" approach.
#'
#' @param reads A \code{\link{GRanges-class}} object with read fragments. NOTE: Complexity estimation relies on duplicate reads and therefore the duplicates have to be present in the input.
#' @return A \code{list} with estimated complexity values and plots.
#' @importFrom stats coefficients nls predict
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

	## Complexity estimation with Michaelis-Menten
	complexity.MM <- NA
	complexity.MM.ggplt <- NA
	vm.init <- quantile(df$y, 1)
	k.init <- quantile(df$x, .25)
	tC <- tryCatch({
		complexity.MM.fit <- stats::nls(y ~ vm * x/(k+x), data=df, start=list(vm=vm.init, k=k.init))
		complexity.MM <- as.numeric(stats::coefficients(complexity.MM.fit)[1])
		max.x <- max(0.9 * stats::coefficients(complexity.MM.fit)['k'] / (1-0.9), max(total.reads))
		x <- seq(from=0, to=max.x, length.out=1000)
		df.fit <- data.frame(x=x, y=stats::predict(complexity.MM.fit, data.frame(x)))
		complexity.MM.ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y'), size=5) + geom_line(data=df.fit, mapping=aes_string(x='x', y='y')) + xlab('total number of reads') + ylab('unique reads') + theme_bw() + ggtitle('Michaelis-Menten complexity estimation')
	}, error = function(err) {
		# warning("Complexity estimation with Michaelis-Menten failed.")
	})

	rl <- list(complexity=c(MM=complexity.MM), MM.ggplt=complexity.MM.ggplt)
	return(rl)
}


