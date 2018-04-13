


#' Filter segments by minimal size
#'
#' \code{filterSegments} filters out segments below a specified minimal segment size. This can be useful to get rid of boundary effects from the Hidden Markov approach.
#'
#' @param segments A \code{\link{GRanges-class}} object.
#' @param min.seg.width The minimum segment width in base-pairs.
#' @return The input \code{model} with adjusted segments.
#' @author Aaron Taudt
#' @export
#'@examples
#'## Load an HMM
#'file <- list.files(system.file("extdata", "primary-lung", "hmms",
#'                   package="AneuFinderData"), full.names=TRUE)
#'hmm <- loadFromFiles(file)[[1]]
#'## Check number of segments before and after filtering
#'length(hmm$segments)
#'hmm$segments <- filterSegments(hmm$segments, min.seg.width=2*width(hmm$bins)[1])
#'length(hmm$segments)
#'
filterSegments <- function(segments, min.seg.width) {
	
	if (is.null(segments)) {
		return(NULL)
	}
	if (min.seg.width<=0) {
		return(segments)
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
	segments.df <- as.data.frame(segments)
	segments.df <- collapseBins(segments.df, column2collapseBy='state', columns2drop=c('width'))
	segments.filtered <- as(segments.df, 'GRanges')
	seqlevels(segments.filtered) <- seqlevels(segments) # correct order after as()
	seqlengths(segments.filtered) <- seqlengths(segments)
	return(segments)
}

#' Get SCE coordinates
#'
#' Extracts the coordinates of a sister chromatid exchanges (SCE) from an \code{\link{aneuBiHMM}} object.
#'
#' @param model An \code{\link{aneuBiHMM}} object.
#' @param resolution An integer vector specifying the resolution at bin level at which to scan for SCE events.
#' @param min.segwidth Segments below this width will be removed before scanning for SCE events.
#' @param fragments A \code{\link{GRanges-class}} object with read fragments or a file that contains such an object. These reads will be used for fine mapping of the SCE events.
#' @return A \code{\link{GRanges-class}} object containing the SCE coordinates.
#' @author Aaron Taudt
#' @export
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the BAM file into bin size 1Mp
#'binned <- binReads(bedfile, assembly='hg19', binsize=1e6,
#'                   chromosomes=c(1:22,'X','Y'), pairedEndReads=TRUE)
#'## Fit the Hidden Markov Model
#'## Find copy-numbers
#'model <- findCNVs.strandseq(binned[[1]])
#'## Find sister chromatid exchanges
#'model$sce <- getSCEcoordinates(model)
#'print(model$sce)
#'plot(model)
#'
getSCEcoordinates <- function(model, resolution=c(3,6), min.segwidth=2, fragments=NULL) {

	if (class(model) != "aneuBiHMM") {
		stop("argument 'model' requires an aneuBiHMM object")
	}
	if (is.null(levels(model$bins$state))) {
		sce <- GRanges()
		return(sce)
	}
	multiplicity <- suppressWarnings( initializeStates(levels(model$bins$state))$multiplicity )

	## Merge '0-somy' and 'zero-inflation'
	bins <- model$bins
	# Add 0-somy and zero-inflation factor levels in case they are not there
	bins$state <- factor(bins$state, levels=unique(c('zero-inflation','0-somy',levels(bins$state))))
	bins$pstate <- factor(bins$pstate, levels=unique(c('zero-inflation','0-somy',levels(bins$pstate))))
	bins$mstate <- factor(bins$mstate, levels=unique(c('zero-inflation','0-somy',levels(bins$mstate))))
	bins$state[bins$state=='zero-inflation'] <- '0-somy'
	bins$mstate[bins$mstate=='zero-inflation'] <- '0-somy'
	bins$pstate[bins$pstate=='zero-inflation'] <- '0-somy'

	## Remove bins with 0-somy in both strands
	bins <- bins[bins$state != 'zero-inflation' & bins$state != '0-somy']
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
	sce <- unlist(sce, use.names=FALSE)
	sce <- sort(sce)
	sce <- reduce(sce)

	return(sce)
}


