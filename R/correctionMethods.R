

#' Mappability correction
#'
#' Correct a list of \code{\link{binned.data}} by mappability.
#'
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @param reference A file or \code{\link{GRanges}} with aligned reads.
#' @param same.binsize If \code{TRUE} the mappability correction will only be calculated once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects describe the same genome at the same binsize.
#' @return A \code{list} with \code{\link{binned.data}} objects with adjusted read counts.
#' @author Aaron Taudt
#' @inheritParams bam2GRanges
#' @inheritParams bed2GRanges
#'
#'
correctMappability <- function(binned.data.list, same.binsize, reference, assembly, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000) {

	binned.data.list <- loadGRangesFromFiles(binned.data.list)
	same.binsize.calculated <- FALSE
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]

		## Calculate GC content per bin
		if (same.binsize & !same.binsize.calculated | !same.binsize) {
			ptm <- startTimedMessage("Calculating mappability per bin ...")
			refbin <- binReads(file=reference, assembly=assembly, chromosomes=seqlevels(binned.data), pairedEndReads=pairedEndReads, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, max.fragment.width=max.fragment.width, binsizes=NULL, reads.per.bin=NULL, bins=list('ref'=binned.data), save.as.RData=FALSE, calc.complexity=FALSE)[[1]]
			## Check if seqlengths of data and mappability correction are consistent
			chromlengths <- seqlengths(binned.data)
			chroms <- names(chromlengths)
			# Compare
			compare <- chromlengths[chroms] == seqlengths(refbin)[chroms]
			if (any(compare==FALSE, na.rm=TRUE)) {
				warning(paste0(attr(binned.data,'ID'),": Chromosome lengths differ between binned data and 'reference'. Mappability correction skipped. Please use the correct genome for option 'reference'."))
				binned.data.list[[i1]] <- binned.data
				next
			}

			## Make the mappability correction vector
			tab <- table(refbin$counts)
			refbin.maxcount <- as.numeric(names(which.max(tab[as.numeric(names(tab))>0])))
			mappability <- refbin$counts / refbin.maxcount
			mappability[mappability==0] <- 1
			
			
			same.binsize.calculated <- TRUE
			stopTimedMessage(ptm)
		}
		binned.data$mappability <- mappability

		### GC correction ###
		ptm <- startTimedMessage("Mappability correction ...")
		counts <- binned.data$counts / binned.data$mappability
		mcounts <- binned.data$mcounts / binned.data$mappability
		pcounts <- binned.data$pcounts / binned.data$mappability
		## Correction factors
		binned.data$counts <- as.integer(round(counts))
		binned.data$mcounts <- as.integer(round(mcounts))
		binned.data$pcounts <- as.integer(round(pcounts))
		binned.data$counts[binned.data$counts<0] <- 0
		binned.data$pcounts[binned.data$pcounts<0] <- 0
		binned.data$mcounts[binned.data$mcounts<0] <- 0

		### Quality measures ###
		## Spikyness
		attr(binned.data, 'spikiness') <- qc.spikiness(binned.data$counts)
		## Shannon entropy
		attr(binned.data, 'shannon.entropy') <- qc.entropy(binned.data$counts)

		binned.data.list[[i1]] <- binned.data
	}
	return(binned.data.list)
}

#' GC correction
#'
#' Correct a list of \code{\link{binned.data}} by GC content.
#'
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
#' @param same.binsize If \code{TRUE} the GC content will only be calculated once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects describe the same genome at the same binsize.
#' @return A \code{list} with \code{\link{binned.data}} objects with adjusted read counts.
#' @author Aaron Taudt
#' @importFrom Biostrings Views alphabetFrequency
#' @importFrom stats lm predict
#' @export
#'@examples
#'## Get a BED file, bin it and run GC correction
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'plot(binned[[1]], type=1)
#'if (require(BSgenome.Mmusculus.UCSC.mm10)) {
#'  binned.GC <- correctGC(list(binned[[1]]), GC.BSgenome=BSgenome.Mmusculus.UCSC.mm10)
#'  plot(binned.GC[[1]], type=1)
#'}
#'
correctGC <- function(binned.data.list, GC.BSgenome, same.binsize=FALSE) {

	binned.data.list <- loadGRangesFromFiles(binned.data.list)
	same.binsize.calculated <- FALSE
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]

		## Check if seqlengths of data and GC.correction are consistent
		# Replace 1->chr1 if necessary
		chromlengths <- seqlengths(binned.data)
		chroms <- names(chromlengths)
		mask <- !grepl('^chr', chroms)
		chroms[mask] <- paste0('chr',chroms[mask])
		names(chromlengths) <- chroms
		# Compare
		compare <- chromlengths[chroms] == seqlengths(GC.BSgenome)[chroms]
		if (any(compare==FALSE, na.rm=TRUE)) {
			warning(paste0(attr(binned.data,'ID'),": Chromosome lengths differ between binned data and 'GC.BSgenome'. GC correction skipped. Please use the correct genome for option 'GC.BSgenome'."))
			binned.data.list[[i1]] <- binned.data
			next
		}

		## Calculate GC content per bin
		if (same.binsize & !same.binsize.calculated | !same.binsize) {
			ptm <- startTimedMessage("Calculating GC content per bin ...")
			GC.content <- list()
			for (chrom in seqlevels(binned.data)) {
				if (!grepl('^chr',chrom)) {
					chr <- paste0('chr',chrom)
				} else {
					chr <- chrom
				}
				if (chr %in% seqlevels(GC.BSgenome)) {
					view <- Biostrings::Views(GC.BSgenome[[chr]], ranges(binned.data)[seqnames(binned.data)==chrom])
					freq <- Biostrings::alphabetFrequency(view, as.prob = TRUE, baseOnly=TRUE)
					if (nrow(freq) > 1) {
						GC.content[[as.character(chrom)]] <- rowSums(freq[, c("G","C")])
					} else {
						GC.content[[as.character(chrom)]] <- sum(freq[, c("G","C")])
					}
				} else {
					GC.content[[as.character(chrom)]] <- rep(NA, length(binned.data[seqnames(binned.data)==chrom]))
					warning(paste0(attr(binned.data,'ID'),": No sequence information for chromosome ",chr," available."))
				}
			}
			GC.content <- unlist(GC.content)
			same.binsize.calculated <- TRUE
			stopTimedMessage(ptm)
		}
		binned.data$GC <- GC.content

		### GC correction ###
		ptm <- startTimedMessage("GC correction ...")
		counts <- binned.data$counts
		mcounts <- binned.data$mcounts
		pcounts <- binned.data$pcounts
		## Correction factors
		gc.categories <- seq(from=0, to=1, length=20)
		intervals.per.bin <- findInterval(binned.data$GC, gc.categories)
		intervals <- sort(unique(intervals.per.bin))
		mean.counts.global <- mean(binned.data$counts, trim=0.05)
		correction.factors <- NULL
		weights <- NULL
		for (interval in intervals) {
			mask <- intervals.per.bin==interval
			counts.with.same.GC <- binned.data$counts[mask]
			weights[as.character(gc.categories[interval])] <- length(counts.with.same.GC)
			mean.counts.with.same.GC <- mean(counts.with.same.GC, na.rm=TRUE, trim=0.05)
			if (mean.counts.with.same.GC == 0) {
				correction.factor <- 0
			} else {
				correction.factor <-  mean.counts.global / mean.counts.with.same.GC
			}
			correction.factors[as.character(gc.categories[interval])] <- correction.factor
		}
		## Fit x^2 to correction.factors
		y <- correction.factors[-1][correction.factors[-1]<10]
		x <- as.numeric(names(y))
		w <- weights[-1][correction.factors[-1]<10]
		df <- data.frame(x,y,weight=w)
		weight <- w	# dummy assignment to pass R CMD check, doesn't affect the fit
		fit <- stats::lm(y ~ poly(x, 2, raw=TRUE), data=df, weights=weight)
		fitted.correction.factors <- stats::predict(fit, data.frame(x=gc.categories[intervals]))
		names(fitted.correction.factors) <- intervals
		for (interval in intervals) {
			mask <- which(intervals.per.bin==interval)
			correction.factor <- fitted.correction.factors[as.character(interval)]
			counts[mask] <- counts[mask] * correction.factor
			mcounts[mask] <- mcounts[mask] * correction.factor
			pcounts[mask] <- pcounts[mask] * correction.factor
		}
		binned.data$counts <- as.integer(round(counts))
		binned.data$mcounts <- as.integer(round(mcounts))
		binned.data$pcounts <- as.integer(round(pcounts))
		binned.data$counts[binned.data$counts<0] <- 0
		binned.data$pcounts[binned.data$pcounts<0] <- 0
		binned.data$mcounts[binned.data$mcounts<0] <- 0

		# Produce fit to check
		ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y', size='weight')) + geom_line(aes_string(x='x', y='y'), data=data.frame(x=gc.categories[intervals], y=fitted.correction.factors)) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('correction factor')
# 		attr(binned.data, 'GC.correction.ggplt') <- ggplt # do not append, ridiculously inflates disk usage
		stopTimedMessage(ptm)

		### Quality measures ###
		## Spikyness
		attr(binned.data, 'spikiness') <- qc.spikiness(binned.data$counts)
		## Shannon entropy
		attr(binned.data, 'shannon.entropy') <- qc.entropy(binned.data$counts)

		binned.data.list[[i1]] <- binned.data
	}
	return(binned.data.list)
}


