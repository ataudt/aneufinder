#' Get IDs of a subset of models
#'
#' Get the IDs of models that have a certain CNV profile. The result will be \code{TRUE} for models that overlap all specified ranges in \code{profile} by at least one base pair with the correct state.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param profile A \code{\link{GRanges-class}} object with metadata column 'expected.state' and optionally columns 'expected.mstate' and 'expected.pstate'.
#' @return A named logical vector with \code{TRUE} for all models that are concordant with the given \code{profile}.
#' @export
#'
#'@examples
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get all files that have a 3-somy on chromosome 1 and 4-somy on chromosome 2
#'profile <- GRanges(seqnames=c('1','2'), ranges=IRanges(start=c(1,1), end=c(195471971,182113224)),
#'                   expected.state=c('3-somy','4-somy'))
#'ids <- subsetByCNVprofile(lung.files, profile)
#'print(which(ids))
#'
subsetByCNVprofile <- function(hmms, profile) {
	
	is.concordant <- logical()
	hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
	for (hmm in hmms) {
		segments <- hmm$segments
		if (!is.null(segments)) {
			mask <- TRUE
			for (iprof in 1:length(profile)) {
				if (!is.null(profile$expected.state)) {
					exp.state <- profile$expected.state[iprof]
					subsegs <- segments[segments$state==exp.state]
					ov <- findOverlaps(subsegs, profile[iprof])
					if (length(ov)==0) { mask <- FALSE }
				}
				if (!is.null(profile$expected.pstate)) {
					exp.state <- profile$expected.pstate[iprof]
					subsegs <- segments[segments$state==exp.state]
					ov <- findOverlaps(subsegs, profile[iprof])
					if (length(ov)==0) { mask <- FALSE }
				}
				if (!is.null(profile$expected.mstate)) {
					exp.state <- profile$expected.mstate[iprof]
					subsegs <- segments[segments$state==exp.state]
					ov <- findOverlaps(subsegs, profile[iprof])
					if (length(ov)==0) { mask <- FALSE }
				}
			}
		} else {
			mask <- FALSE
		}
		is.concordant[hmm$ID] <- mask
	}
	return(is.concordant)

}


# #' Concatenate GRanges in range
# #'
# #' Concatenate \code{\link{GRanges-class}} objects in range.
# #'
# #' @param grlist A list with \code{\link{GRanges-class}} objects or files that contain such objects.
# #' @param chromosome,start,end Coordinates of the range.
# #' @export
# concGRangesInRange <- function(grlist, chromosome, start, end) {
# 
# 	range <- GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end))
# 	
# 	concGR <- GRangesList()
# 	for (i1 in 1:length(grlist)) {
# 		gr <- suppressMessages( loadFromFiles(grlist[[i1]], check.class='GRanges')[[1]] )
# 		concGR[[i1]] <- subsetByOverlaps(gr, range)
# 	}
# 	concGR <- unlist(concGR, use.names=FALSE)
# 	return(concGR)
# 
# }
# 
# 
# #' Smooth binning of read fragments
# #'
# #' The function bins read fragments with a fixed-length sliding window.
# #'
# #' @param gr A \code{\link{GRanges-class}} object. The GRanges is assumed to have only data for one chromosome.
# #' @param binsize The size in base-pairs of the fixed-length window.
# #' @param stepsize The number of base-pairs for sliding the window.
# #' @param start,end Coordinates of the binning region. If \code{NULL} coordinates will be determined from the data.
# smoothBinning <- function(gr, binsize, stepsize=binsize/10, start=NULL, end=NULL) {
# 
# 	if (length(unique(seqnames(gr))) > 1) {
# 		stop("argument 'gr' must contain data for only one chromosome")
# 	}
# 	if (is.null(start)) {
# 		start <- start(gr)[1]
# 	}
# 	if (is.null(end)) {
# 		end <- end(gr)[length(gr)]
# 	}
# 
# 	offsets <- seq(from=0, to=binsize-1, by=stepsize)
# 	binned.list <- GRangesList()
# 	for (ioff in 1:length(offsets)) {
# 		## Make the GRanges for the binned read counts ##
# 		startpos <- start + offsets[ioff]
# 		endpos <- end
# 		starts <- seq(from=startpos, to=endpos, by=binsize)
# 		ends <- starts[-1]-1
# 		starts <- starts[-length(starts)]
# 		binned <- GRanges(seqnames=unique(seqnames(gr)), ranges=IRanges(start=starts, end=ends))
# 		## Count overlaps
# 		binned$counts <- countOverlaps(binned,gr)
# 		binned.list[[ioff]] <- binned
# 	}
# 	binned <- sort(unlist(binned.list, use.names=FALSE))
# 	return(binned)
# 	
# }
# 
