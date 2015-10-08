#' Get IDs of a subset of models
#'
#' Get the IDs of models that have a certain CNV profile.
#'
#' @param hmm.list A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param profile A \code{\link{GRanges}} object with metadata column 'expected.state' and optionally columns 'expected.mstate' and 'expected.pstate'.
#' @return A vector with the IDs of the models that are concordant with the given \code{profile}.
#' @export
subsetByCNVprofile <- function(hmm.list, profile) {
	
	hmm.list <- loadHmmsFromFiles(hmm.list)
	is.concordant <- sapply(hmm.list, function(hmm) {
		ind <- findOverlaps(profile, hmm$segments, select='first')
		mask <- TRUE
		if (!is.null(profile$expected.state)) {
			mask <- mask & all(hmm$segments$state[ind] == profile$expected.state)
		}
		if (!is.null(profile$expected.mstate)) {
			mask <- mask & all(hmm$segments$mstate[ind] == profile$expected.mstate)
		}
		if (!is.null(profile$expected.pstate)) {
			mask <- mask & all(hmm$segments$pstate[ind] == profile$expected.pstate)
		}
		return(mask)
	})
	ids <- sapply(hmm.list, '[[', 'ID')
	con.ids <- ids[is.concordant]
	return(con.ids)
}


#' Concatenate GRanges in range
#'
#' Concatenate \code{\link{GRanges}} objects in range.
#'
#' @param grlist A list with \code{\link{GRanges}} objects or files that contain such objects.
#' @param chromosome,start,end Coordinates of the range.
#' @export
concGRangesInRange <- function(grlist, chromosome, start, end) {

	range <- GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end))
	
	concGR <- GRangesList()
	for (i1 in 1:length(grlist)) {
		gr <- suppressMessages( loadGRangesFromFiles(grlist[[i1]])[[1]] )
		concGR[[i1]] <- subsetByOverlaps(gr, range)
	}
	concGR <- unlist(concGR)
	return(concGR)

}


#' Smooth binning of read fragments
#'
#' The function bins read fragments with a fixed-length sliding window.
#'
#' @param gr A \code{\link{GRanges}} object. The GRanges is assumed to have only data for one chromosome.
#' @param binsize The size in base-pairs of the fixed-length window.
#' @param stepsize The number of base-pairs for sliding the window.
#' @param start,end Coordinates of the binning region. If \code{NULL} coordinates will be determined from the data.
#' @export
smoothBinning <- function(gr, binsize, stepsize=binsize/10, start=NULL, end=NULL) {

	if (length(unique(seqnames(gr))) > 1) {
		stop("argument 'gr' must contain data for only one chromosome")
	}
	if (is.null(start)) {
		start <- start(gr)[1]
	}
	if (is.null(end)) {
		end <- end(gr)[length(gr)]
	}

	offsets <- seq(from=0, to=binsize-1, by=stepsize)
	binned.list <- GRangesList()
	for (ioff in 1:length(offsets)) {
		## Make the GRanges for the binned read counts ##
		startpos <- start + offsets[ioff]
		endpos <- end
		starts <- seq(from=startpos, to=endpos, by=binsize)
		ends <- starts[-1]-1
		starts <- starts[-length(starts)]
		binned <- GRanges(seqnames=unique(seqnames(gr)), ranges=IRanges(start=starts, end=ends))
		## Count overlaps
		binned$counts <- countOverlaps(binned,gr)
		binned.list[[ioff]] <- binned
	}
	binned <- sort(unlist(binned.list))
	return(binned)
	
}

