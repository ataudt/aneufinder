#' Calculate deltaWs
#'
#' This function will calculate deltaWs from a \code{\link{GRanges}} object with read fragments.
#'
#' @param frags A \code{\link{GRanges}} with read fragments (see \code{\link{bam2GRanges}}).
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

	return(frags.new)

}
