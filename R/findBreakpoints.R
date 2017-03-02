#' Find breakpoints
#' 
#' Breakpoint detection is done via a dynamic windowing approach on read resolution.
#' 
#' @param model An \code{\link{aneuBiHMM}} object or a file that contains such an object.
#' @param fragments A \code{\link[GenomicRanges]{GRanges}} object with read fragments.
#' @param breakpoint.quantile A quantile cutoff between 0 and 1 for breakpoint detection. Higher values will result in higher precision but lower sensitivity.
#' @author Aaron Taudt, David Porubsky, Ashley Sanders
findBreakpoints <- function(model, fragments, breakpoint.quantile=0.99) {
  
  ## Load model
  model <- loadFromFiles(model, check.class = 'aneuBiHMM')[[1]]
  if (is.character(fragments)) {
    if (!file.exists(fragments)) {
      stop("Could not find file ", fragments)
    }
  }
  
  ## Window size for deltaWCalculator
  reads.per.window <- as.integer(mean(model$bins$counts))
  
  ## Raw breakpoints
  ptm <- startTimedMessage("Calculating deltaWs ...")
	fragments <- suppressWarnings( deltaWCalculator(fragments, reads.per.window=reads.per.window) )
	fragments$cdf <- ecdf(fragments$deltaW)(fragments$deltaW)
	stopTimedMessage(ptm)

	## Make peak list
  ptm <- startTimedMessage("Getting peak numbers ...")
	mask <- fragments$cdf >= breakpoint.quantile
	rlemask <- rle(mask)
	rlegroup <- rlemask
	rlemaskT <- rlemask$values==TRUE
	rlegroup$values[rlemaskT] <- cumsum(rlemask$values[rlemaskT])
	fragments$peakGroup <- inverse.rle(rlegroup)
	stopTimedMessage(ptm)
	
	## Select the maximum deltaW read(s) within each peak
  ptm <- startTimedMessage("Peak summit ...")
	peaks.unrefined <- fragments[fragments$peakGroup>0]
	peaksummits <- GRangesList()
	for (group in unique(peaks.unrefined$peakGroup)) {
	  peak <- peaks.unrefined[peaks.unrefined$peakGroup == group]
    maxpeakdeltaW <- max(peak$deltaW)
    peaksummit <- peak[peak$deltaW == max(peak$deltaW)]
    strand(peaksummit) <- '*'
    rpeaksummit <- range(peaksummit)
    mcols(rpeaksummit)[c('deltaW', 'peakGroup', 'cdf')] <- mcols(peaksummit)[c('deltaW', 'peakGroup', 'cdf')][1,]
    rpeaksummit$numReadsInPeak <- length(peak)
    peaksummits[[as.character(group)]] <- rpeaksummit
	}
	peaks <- unlist(peaksummits, use.names = FALSE)
	stopTimedMessage(ptm)
	
	## Genotyping for each strand
  ptm <- startTimedMessage("Genotyping ...")
	# Interval between breaks/peaksummits
	intervals <- gaps(peaks)
	intervals <- intervals[strand(intervals) == '*']
	# Counts per interval
	strand <- '-'
	strand(intervals) <- strand
	intervals$mcounts <- countOverlaps(intervals, fragments)
	strand <- '+'
	strand(intervals) <- strand
	intervals$pcounts <- countOverlaps(intervals, fragments)
	# Normalize
	binsize <- mean(width(model$bins))
	intervals$mcounts.n <- intervals$mcounts / width(intervals) * binsize
	intervals$pcounts.n <- intervals$pcounts / width(intervals) * binsize
	intervals$mCN <- round(intervals$mcounts.n / model$distributions$minus['1-somy', 'mu'])
	intervals$pCN <- round(intervals$pcounts.n / model$distributions$plus['1-somy', 'mu'])
	intervals$CN <- paste(intervals$mCN, intervals$pCN)
	# Split by chromosome and remove consecutive intervals with the same copy number
	peaks$CN.from <- NA
	peaks$CN.to <- NA
	peaks.list <- GRangesList()
	for (chrom in seqlevels(peaks)) {
	  peaks.chrom <- peaks[seqnames(peaks) == chrom]
	  intervals.chrom <- intervals[seqnames(intervals) == chrom]
	  peaks.chrom$CN.from <- intervals.chrom$CN[-length(intervals.chrom)]
	  peaks.chrom$CN.to <- intervals.chrom$CN[-1]
	  peaks.list[[chrom]] <- peaks.chrom
	}
	peaks <- unlist(peaks.list, use.names = FALSE)
	peaks <- peaks[peaks$CN.from != peaks$CN.to]
	stopTimedMessage(ptm)
	
	return(peaks)
	
	
}



#' Calculate deltaWs
#'
#' This function will calculate deltaWs from a \code{\link[GenomicRanges]{GRanges}} object with read fragments.
#'
#' @param frags A \code{\link[GenomicRanges]{GRanges}} with read fragments (see \code{\link{bam2GRanges}}).
#' @param reads.per.window Number of reads in each dynamic window.
#' @return The input \code{frags} with additional meta-data columns.
#' @import GenomicRanges
#' @importFrom BiocGenerics as.vector
#' @author Aaron Taudt, David Porubsky, Ashley Sanders
#' @export
deltaWCalculator <- function(frags, reads.per.window=10) {

	if (reads.per.window == 0) {
		stop("'reads.per.window' must be >= 1")
	}
	if (reads.per.window < 10) {
		warning("'reads.per.window' should at least be 10")
	}
	if (is.character(frags)) {
		frags <- loadFromFiles(frags, check.class='GRanges')[[1]]
	}
	frags.split <- split(frags, seqnames(frags))
	reads.per.chrom <- sapply(frags.split, length)
	chroms2parse <- names(reads.per.chrom)[reads.per.chrom>2*reads.per.window]
	chroms2skip <- setdiff(names(reads.per.chrom),chroms2parse)
	if (length(chroms2skip)>0) {
		warning(paste0("Not parsing chromosomes ",paste(chroms2skip, collapse=',')," because they do not have enough reads."))
	}
	if (length(chroms2parse)==0) {
		warning("None of the specified chromosomes has enough reads. Doing nothing.")
		return(GRanges())
	}

	frags.new <- GRangesList()
	for (chrom in chroms2parse) {
		f <- frags.split[[chrom]]
		f <- f[order(start(f))]
		f$pcsum <- cumsum(strand(f)=='+')
		f$mcsum <- cumsum(strand(f)=='-')
		f$preads <- c(rep(NA,reads.per.window),diff(BiocGenerics::as.vector(f$pcsum),lag=reads.per.window))
		f$mreads <- c(rep(NA,reads.per.window),diff(BiocGenerics::as.vector(f$mcsum),lag=reads.per.window))
		f$deltaW <- abs(c(diff(f$preads,lag=reads.per.window),rep(NA,reads.per.window)))
		# Shift deltaWs to region between reads
		start.f <- end(f)
		end.f <- c(start(f)[-1],seqlengths(frags)[chrom])
		mask <- start.f < end.f
		f <- f[mask]
		start(f) <- start.f[mask]
		end(f) <- end.f[mask]
		frags.new[[chrom]] <- f
	}
	frags.new <- unlist(frags.new)
	names(frags.new) <- NULL
	# Replace NAs with 0 to avoid problems in downstream functions
	frags.new$deltaW[is.na(frags.new$deltaW)] <- 0
	frags.new$mreads[is.na(frags.new$mreads)] <- 0
	frags.new$preads[is.na(frags.new$preads)] <- 0
	## Remove unneeded columns
	frags.new$pcsum <- NULL
	frags.new$mcsum <- NULL

	return(frags.new)

}
