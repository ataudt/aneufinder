

#' Find sister chromatid exchanges
#'
#' \code{findSCEs} classifies the binned read counts into several states which represent the number of chromatids on each strand.
#'
#' \code{findSCEs} uses a Hidden Markov Model to classify the binned read counts: state 'zero-inflation' with a delta function as emission densitiy (only zero read counts), 'nullsomy' with geometric distribution, 'monosomy','disomy','trisomy','tetrasomy' and 'multisomy' with negative binomials (see \code{\link{dnbinom}}) as emission densities. A expectation-maximization (EM) algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method. TODO: insert paper
#' @author Aaron Taudt
#' @inheritParams univariate.findCNVs
#' @inheritParams bivariate.findCNVs
#' @inheritParams getSCEcoordinates
#' @return An \code{\link{aneuBiHMM}} object.
#' @examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata/BB140820_I_002.bam", package="aneufinder")
#'## Bin the BAM file into bin size 200000bp
#'binned.data <- bam2binned(bamfile, binsize=200000, chromosomes=c(1:22,'X','Y'), save.as.RData=FALSE)
#'## Fit the Hidden Markov Model
#'model <- findSCEs(binned.data, ID=basename(bamfile), eps=0.1, max.time=60)
#'## Check the fit
#'plot(model, type='histogram')
#' @export
findSCEs <- function(binned.data, ID=NULL, eps=0.1, init="standard", max.time=-1, max.iter=1000, num.trials=5, eps.try=10*eps, num.threads=1, count.cutoff.quantile=0.999, strand='*', allow.odd.states=TRUE, states=c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"), most.frequent.state="monosomy", algorithm="EM", initial.params=NULL) {

	## Intercept user input
	if (class(binned.data) != 'GRanges') {
		binned.data <- get(load(binned.data))
		if (class(binned.data) != 'GRanges') stop("argument 'binned.data' expects a GRanges with meta-column 'counts' or a file that contains such an object")
	}
	if (is.null(ID)) {
		ID <- attr(binned.data, 'ID')
	}

	## Print some stuff
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	message("Find CNVs for ID = ",ID, ":")

	model <- bivariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, count.cutoff.quantile=count.cutoff.quantile, allow.odd.states=allow.odd.states, states=states, most.frequent.state=most.frequent.state, algorithm=algorithm, initial.params=initial.params)
	
# 	## Find CNV calls for offset counts using the parameters from the normal run
# 	offsets <- setdiff(names(attr(binned.data,'offset.counts')), 0)
# 	if (!is.null(offsets)) {
# 		offset.models <- list()
# 		for (ioff in offsets) {
# 			message(paste0("Finding SCE for offset ",ioff))
# 			off.counts <- attr(binned.data,'offset.counts')[[as.character(ioff)]]
# 			off.binned.data <- binned.data
# 			mcols(off.binned.data)[names(mcols(binned.data)) %in% names(off.counts)] <- as(off.counts, 'DataFrame')
# 			off.model <- suppressMessages( bivariate.findCNVs(off.binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, count.cutoff.quantile=count.cutoff.quantile, allow.odd.states=allow.odd.states, states=states, most.frequent.state=most.frequent.state, algorithm='baumWelch', initial.params=model) )
# 			offset.models[[as.character(ioff)]] <- off.model
# 		}
# 	}
# 	return(offset.models)

	attr(model, 'call') <- call
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(model)

}


#' Filter segments by minimal size
#'
#' \code{filterSegments} filters out segments below a specified minimal segment size. This can be useful to get rid of boundary effects from the Hidden Markov approach.
#'
#' @param model A \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} object.
#' @param min.seg.width.bp The minimum segment width in base-pairs.
#' @param min.seg.width.binsize The minimum segment width in bins.
#' @return The input \code{model} with adjusted segments.
#' @author Aaron Taudt
#' @export
filterSegments <- function(model, min.seg.width.binsize=2, min.seg.width.bp=NULL) {
	
	if (!is.null(min.seg.width.bp)) {
		min.seg.width <- min.seg.width.bp
	} else {
		min.seg.width <- min.seg.width.binsize*width(model$bins)[1]
	}
	if (min.seg.width<=0) {
		return(model)
	}

	segments <- model$segments
	if (is.null(segments)) {
		return(model)
	}
	replace.index <- which(width(segments) < min.seg.width)
	repl.segments <- segments[replace.index]
	keep.index <- which(width(segments) >= min.seg.width)
	keep.segments <- segments[keep.index]
	nearest.index <- nearest(repl.segments, keep.segments)
	na.mask <- is.na(nearest.index)
	nearest.index <- nearest.index[!na.mask]
	replace.index <- replace.index[!na.mask]
	if (length(nearest.index)>0) {
		nearest.keep.segments <- keep.segments[nearest.index]
		segments$state[replace.index] <- nearest.keep.segments$state
		segments$mstate[replace.index] <- nearest.keep.segments$mstate
		segments$pstate[replace.index] <- nearest.keep.segments$pstate
	}
	model$segments <- segments

	return(model)
}

#' Get SCE coordinates
#'
#' Extracts the coordinates of a sister chromatid exchanges (SCE) from an \code{\link{aneuBiHMM}} object.
#'
#' @param model An \code{\link{aneuBiHMM}} object.
#' @param resolution An integer vector specifying the resolution at bin level at which to scan for SCE events.
#' @param min.segwidth Minimum segment length in bins when scanning for SCE events.
#' @param fragments A \code{\link{GRanges}} object with read fragments or a file that contains such an object. These reads will be used for fine mapping of the SCE events.
#' @param min.reads Minimum number of reads required for refinement.
#' @return A \code{\link{GRanges}} object containing the SCE coordinates.
#' @author Aaron Taudt
#' @export
getSCEcoordinates <- function(model, resolution=c(3,6), min.segwidth=2, fragments=NULL, min.reads=50) {

	if (is.null(levels(model$bins$state))) {
		sce <- GRanges()
		return(sce)
	}
	multiplicity <- initializeStates(levels(model$bins$state))$multiplicity

	## Merge 'nullsomy' and 'zero-inflation'
	bins <- model$bins
	bins$state[bins$state=='zero-inflation'] <- 'nullsomy'
	bins$mstate[bins$mstate=='zero-inflation'] <- 'nullsomy'
	bins$pstate[bins$pstate=='zero-inflation'] <- 'nullsomy'

	## Remove bins with nullsomy in both strands
	bins <- bins[bins$state != 'zero-inflation' & bins$state != 'nullsomy']
	## Remove bins below minimum segment size
	minsegs <- model$segments[width(model$segments) >= min.segwidth*width(model$bins)[1]]
	bins <- subsetByOverlaps(bins, minsegs)

	## Find SCE coordinates for each chromosome
	bins.split <- split(bins, seqnames(bins))
	sce <- GRangesList()
	for (chrom in names(bins.split)) {
		bins.chrom <- bins.split[[chrom]]
		if (length(bins.chrom)>1) {
			multiplicity.minus <- multiplicity[bins.chrom$mstate]
			multiplicity.plus <- multiplicity[bins.chrom$pstate]
			multiplicity.both <- multiplicity[bins.chrom$state]
			for (ires in resolution) {
				diff.m <- c(rep(0,ires),diff(multiplicity.minus,lag=ires))
				diff.p <- c(rep(0,ires),diff(multiplicity.plus,lag=ires))
				index <- which((diff.m > 0 & diff.p < 0) | (diff.m < 0 & diff.p > 0))
				if (length(index)>0) {
					sce.new <- bins.chrom[index]
					start(sce.new) <- start(bins.chrom[index-ires])
					mcols(sce.new) <- NULL
					if (!is.null(sce[[as.character(chrom)]])) {
						sce.new <- setdiff(sce.new, subsetByOverlaps(sce.new, sce[[as.character(chrom)]]))
						sce[[as.character(chrom)]] <- c(sce[[as.character(chrom)]], sce.new)
					} else {
						sce[[as.character(chrom)]] <- sce.new
					}
				}
			}
		}
	}
	sce <- unlist(sce)
	names(sce) <- NULL
	sce <- sort(sce)
	sce <- reduce(sce)

	### Fine mapping of each SCE ###
	if (!is.null(fragments) & length(sce)>0) {
		fragments <- loadGRangesFromFiles(fragments)[[1]]
		starts <- start(sce)
		ends <- end(sce)
		for (isce in 1:length(sce)) {
			ifrags <- subsetByOverlaps(fragments, sce[isce])
			if (length(ifrags)>5*min.reads) {
				f <- ifrags[order(start(ifrags))]
				f$pcsum <- cumsum(strand(f)=='+')
				f$mcsum <- cumsum(strand(f)=='-')
				f$preads <- c(rep(NA,min.reads),diff(BiocGenerics::as.vector(f$pcsum),lag=min.reads))
				f$mreads <- c(rep(NA,min.reads),diff(BiocGenerics::as.vector(f$mcsum),lag=min.reads))
				f$pdelta <- abs(c(diff(f$preads,lag=min.reads),rep(NA,min.reads)))
				ifrags <- f
				indmax <- which(ifrags$pdelta==max(ifrags$pdelta,na.rm=T))
				starts[isce] <- start(ifrags)[indmax[1]]
				ends[isce] <- end(ifrags)[indmax[length(indmax)]]
			}
		}
		sce.fine <- sce
		start(sce.fine) <- starts
		end(sce.fine) <- ends
		return(sce.fine)
	}
	return(sce)
}


