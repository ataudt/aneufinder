#' Make consensus segments
#'
#' Make consensus segments from a list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects.
#'
#' The function will produce a \code{\link[GenomicRanges]{GRanges-class}} object using the \code{GenomicRanges::disjoin} function on all extracted \code{$segment} entries.
#'
#' @param hmms A list of \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects or a character vector of files that contains such objects.
#' @return A \code{\link[GenomicRanges]{GRanges-class}}.
#' @export
#' @examples 
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get consensus segments and states
#'consensusSegments(lung.files)
#' 
consensusSegments <- function(hmms) {

  	## Load the files
  	hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
  
    ## Get segments from list
    segs.list <- GRangesList()
    for (hmm in hmms) {
        if (!is.null(hmm$segments)) {
            segs.list[[hmm$ID]] <- hmm$segments
        }
    }
    ## Consensus template
    consensus <- disjoin(unlist(segs.list, use.names=FALSE))
    constates <- array(NA, dim=c(length(consensus), length(hmms)), dimnames=list(chromosome=as.character(seqnames(consensus)), sample=names(segs.list)))
    for (i1 in 1:length(segs.list)) {
        segs <- segs.list[[i1]]
        splt <- split(segs, mcols(segs)$copy.number)
        mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
        constates[,i1] <- as.numeric(names(splt)[mind])
    }
    consensus$copy.number <- constates
  
  	return(consensus)
}

