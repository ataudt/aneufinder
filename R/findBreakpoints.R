#' Extract breakpoints
#' 
#' Extract breakpoints with confidence intervals from an \code{\link{aneuBiHMM}} object.
#' 
#' Confidence intervals for breakpoints are estimated by going outwards from the breakpoint read by read, and performing a binomial test of getting the observed or a more extreme outcome, given that the reads within the confidence interval belong to the other side of the breakpoint.
#' 
#' @param model An \code{\link{aneuBiHMM}} object or a file that contains such an object.
#' @param fragments A \code{\link{GRanges}} object with read fragments.
#' @param conf Desired confidence interval.
#' @return A \code{\link{GRanges}} with breakpoint coordinates and confidence interals if \code{fragments} was specified.
#' @examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the data into bin size 1Mp
#'readfragments <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'), return.reads=TRUE)
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'## Fit the Hidden Markov Model
#'model <- findCNVs.strandseq(binned[[1]], eps=0.1, max.time=60)
#'## Add confidence intervals
#'breaks <- getBreakpoints(model, readfragments)
#' 
getBreakpoints <- function(model, fragments=NULL, conf=0.99) {
  
    model <- loadFromFiles(model, check.class = c("aneuHMM", "aneuBiHMM"))[[1]]
    
    ## Get breakpoints
    breaks <- model$segments
    mcols(breaks) <- NULL
    end(breaks) <- end(breaks) - 1
    seqlengths(breaks) <- NA # to not have the chromosome end as gap
    breaks <- gaps(breaks)
    breaks <- breaks[strand(breaks)=='*']
    seqlengths(breaks) <- seqlengths(model$segments)
    # Add genotype transition
    ind <- findOverlaps(breaks, model$segments)
    breaks$mstate.left <- factor(NA, levels=levels(model$segments$mstate))
    breaks$mstate.left[ind@from] <- model$segments$mstate[ind@to]
    breaks$pstate.left <- factor(NA, levels=levels(model$segments$pstate))
    breaks$pstate.left[ind@from] <- model$segments$pstate[ind@to]
    breaks$mstate.right <- factor(NA, levels=levels(model$segments$mstate))
    breaks$mstate.right[ind@from] <- model$segments$mstate[ind@to+1]
    breaks$pstate.right <- factor(NA, levels=levels(model$segments$pstate))
    breaks$pstate.right[ind@from] <- model$segments$pstate[ind@to+1]
    
    if (is.null(fragments)) {
        return(breaks)
    }
    
    ## Sort fragments
    seqlevels(fragments) <- seqlevels(model$segments)
    fragments <- sort(fragments, ignore.strand=TRUE)
    fragments$p <- 0
    
    ## Distributions
    if (class(model) == 'aneuHMM') {
        distr <- model$distributions
    } else if (class(model) == 'aneuBiHMM') {
        distr <- model$univariateParams$distributions
    }
    
    ## Do chromosomes one by one
    breaks.conf <- GenomicRanges::GRangesList()
    seqlevels(breaks.conf) <- seqlevels(breaks)
    for (chrom in unique(seqnames(breaks))) {
        message("chrom = ", chrom)
        cbreaks <- breaks[seqnames(breaks) == chrom]
        cfrags <- fragments[seqnames(fragments) == chrom]
        if (length(cbreaks) > 0) {
            for (ibreak in 1:length(cbreaks)) {
                # message(ibreak)
                states <- mcols(cbreaks[ibreak])
                probs <- array(NA, dim=c(2,2), dimnames=list(strand=c('-','+'), direction=c('left','right')))
                probs['-','left'] <- distr[as.character(states[,'mstate.right']),'mu']
                probs['+','left'] <- distr[as.character(states[,'pstate.right']),'mu'] # we need right(!) side probabilities for the left side
                probs['-','right'] <- distr[as.character(states[,'mstate.left']),'mu']
                probs['+','right'] <- distr[as.character(states[,'pstate.left']),'mu'] # we need right(!) side probabilities for the right side
                probs <- sweep(probs, MARGIN = 2, STATS = colSums(probs), FUN = '/')
                # Left side
                ind <- which(start(cfrags) < start(cbreaks)[ibreak])
                if (length(ind) > 0) {
                    ind <- ind[length(ind)]
                    p <- 1
                    numReads <- c('-'=0, '+'=0)
                    if (!any(is.na(probs))) {
                        i1 <- -1
                        while (p > 1-conf) {
                            i1 <- i1+1
                            if (ind-i1 <= 0) {
                                i1 = i1 - 1 
                                break
                            }
                            strandofread <- as.character(strand(cfrags)[ind-i1])
                            strandtocompare <- names(which(probs[, 'left'] >= probs[, 'right']))[1] # directionality of test
                            numReads[strandofread] <- numReads[strandofread] + 1
                            p <- pbinom(q = numReads[strandtocompare], size = sum(numReads), prob = probs[strandofread,'left'])
                        }
                        start(cbreaks)[ibreak] <- start(cfrags)[ind-i1]
                    }
                }
                # Right side
                ind <- which(end(cfrags) > end(cbreaks)[ibreak])
                if (length(ind) > 0) {
                    ind <- ind[1]
                    p <- 1
                    numReads <- c('-'=0, '+'=0)
                    if (!any(is.na(probs))) {
                        i1 <- -1
                        while (p > 1-conf) {
                            i1 <- i1+1
                            if (ind+i1 > length(cfrags)) {
                                i1 = i1 - 1
                                break
                            }
                            strandofread <- as.character(strand(cfrags)[ind+i1])
                            strandtocompare <- names(which(probs[, 'right'] >= probs[, 'left']))[1] # directionality of test
                            numReads[strandofread] <- numReads[strandofread] + 1
                            p <- pbinom(q = numReads[strandtocompare], size = sum(numReads), prob = probs[strandofread,'right'])
                        }
                        end(cbreaks)[ibreak] <- end(cfrags)[ind+i1]
                    }
                }
            }
        }
        breaks.conf[[chrom]] <- cbreaks
    }
    breaks.conf <- unlist(breaks.conf, use.names = FALSE)
    
    # mcols(breaks) <- NULL
    breaks$start.conf <- start(breaks.conf)
    breaks$end.conf <- end(breaks.conf)
    return(breaks)
  
}