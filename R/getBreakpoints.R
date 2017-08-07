#' Extract breakpoints
#' 
#' Extract breakpoints with confidence intervals from an \code{\link{aneuBiHMM}} object.
#' 
#' Confidence intervals for breakpoints are estimated by going outwards from the breakpoint read by read, and performing a test of getting the observed or a more extreme outcome, given that the reads within the confidence interval belong to the other side of the breakpoint.
#' 
#' @param model An \code{\link{aneuBiHMM}} object or a file that contains such an object.
#' @param fragments A \code{\link{GRanges}} object with read fragments or a file that contains such an object.
#' @param confint Desired confidence interval for breakpoints.
#' @return A \code{\link{GRanges}} with breakpoint coordinates and confidence interals if \code{fragments} was specified.
#' @importFrom stats pnbinom pbinom pgeom
#' @export
#' 
#' @examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the data into bin size 1Mp
#'readfragments <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'), reads.return=TRUE)
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'## Fit the Hidden Markov Model
#'model <- findCNVs.strandseq(binned[[1]], eps=0.1, max.time=60)
#'## Add confidence intervals
#'breaks <- getBreakpoints(model, readfragments)
#' 
getBreakpoints <- function(model, fragments=NULL, confint=0.99) {
  
    model <- loadFromFiles(model, check.class = c("aneuHMM", "aneuBiHMM"))[[1]]
    fragments <- loadFromFiles(fragments, check.class = 'GRanges')[[1]]
    binsize <- mean(width(model$bincounts[[1]]))
    
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
    distr <- model$distributions
    if (class(distr) == 'list') {
        distr <- distr[[1]]
    }
    
    breaks.conf <- confidenceIntervals(breaks = breaks, fragments = fragments, distr = distr, confint = confint, binsize = binsize)
    
    breaks$start.conf <- start(breaks.conf)
    breaks$end.conf <- end(breaks.conf)
    return(breaks)
  
}


    
confidenceIntervals <- function(breaks, fragments, distr, confint, binsize) {
  
    ptm <- startTimedMessage("Finding confidence intervals ...")
    min.i1 <- 10
    ## Do chromosomes one by one
    breaks.conf <- GenomicRanges::GRangesList()
    seqlevels(breaks.conf) <- seqlevels(breaks)
    for (chrom in unique(seqnames(breaks))) {
        # message("chrom = ", chrom)
        cbreaks <- breaks[seqnames(breaks) == chrom]
        cbreaks.mod <- cbreaks
        cfrags <- fragments[seqnames(fragments) == chrom]
        if (length(cbreaks) > 0) {
            for (ibreak in 1:length(cbreaks)) {
                # message(ibreak)
              
                states <- as.data.frame(mcols(cbreaks[ibreak]))
                mus <- array(NA, dim=c(2,2), dimnames=list(strand=c('-','+'), direction=c('left','right')))
                mus['-','right'] <- distr[as.character(states[,'mstate.right']),'mu']
                mus['+','right'] <- distr[as.character(states[,'pstate.right']),'mu']
                mus['-','left'] <- distr[as.character(states[,'mstate.left']),'mu']
                mus['+','left'] <- distr[as.character(states[,'pstate.left']),'mu']
                # Directionality of test
                left.is.bigger <- mus[,'left'] >= mus[,'right']
                right.is.bigger <- !left.is.bigger
                
                ## Left side
                ind <- which(start(cfrags) < start(cbreaks)[ibreak])
                if (ibreak > 1) {
                    max.bp <- start(cbreaks[ibreak]) - start(cbreaks[ibreak-1])
                } else {
                    max.bp <- start(cbreaks[ibreak])
                }
                if (length(ind) > 0) {
                    ind <- ind[length(ind)]
                    p <- 1
                    ps <- numeric()
                    ps['0'] <- p
                    numReads <- c('-'=0, '+'=0)
                    if (!any(is.na(mus))) {
                        i1 <- -1
                        ## Helpers for speed improvement
                        start.cfrags <- start(cfrags)
                        start.cbreaks <- start(cbreaks)
                        strand.cfrags <- as.character(strand(cfrags))
                        while (p > 1-confint | i1 <= min.i1) {
                            i1 <- i1+1
                            if (ind-i1 <= 0) {
                                i1 = i1 - 1
                                break
                            }
                            i1.bp <- - (start.cfrags[ind-i1] - start.cbreaks[ibreak])
                            if (i1.bp > max.bp) {
                                i1 = i1 - 1
                                break
                            }
                            strandofread <- strand.cfrags[ind-i1]
                            numReads[strandofread] <- numReads[strandofread] + 1
                            
                            # Minus strand
                            dtype <- distr[as.character(states[,'mstate.right']), 'type']
                            if (dtype == 'dnbinom') {
                                p.minus <- stats::pnbinom(q = numReads['-'] - !right.is.bigger['-'], size = distr[as.character(states[,'mstate.right']), 'size'] * i1.bp/binsize, prob = distr[as.character(states[,'mstate.right']), 'prob'], lower.tail = right.is.bigger['-'])
                            } else if (dtype == 'dgeom') {
                                p.minus <- stats::pgeom(q = numReads['-'] - !right.is.bigger['-'], prob = dgeom.prob(distr[as.character(states[,'mstate.right']), 'mu'] * i1.bp/binsize), lower.tail = right.is.bigger['-'])
                            } else if (dtype == 'delta') {
                                p.minus <- c('-'=as.numeric(numReads['-'] == 0))
                            } else if (dtype == 'dbinom') {
                                p.minus <- stats::pbinom(q = numReads['-'] - !right.is.bigger['-'], size = round(distr[as.character(states[,'mstate.right']), 'size'] * i1.bp/binsize), prob = distr[as.character(states[,'mstate.right']), 'prob'], lower.tail = right.is.bigger['-'])
                            }
                            # Plus strand
                            dtype <- distr[as.character(states[,'pstate.right']), 'type']
                            if (dtype == 'dnbinom') {
                                p.plus <- stats::pnbinom(q = numReads['+'] - !right.is.bigger['+'], size = distr[as.character(states[,'pstate.right']), 'size'] * i1.bp/binsize, prob = distr[as.character(states[,'pstate.right']), 'prob'], lower.tail = right.is.bigger['+'])
                            } else if (dtype == 'dgeom') {
                                p.plus <- stats::pgeom(q = numReads['+'] - !right.is.bigger['+'], prob = dgeom.prob(distr[as.character(states[,'pstate.right']), 'mu'] * i1.bp/binsize), lower.tail = right.is.bigger['+'])
                            } else if (dtype == 'delta') {
                                p.plus <- c('+'=as.numeric(numReads['+'] == 0))
                            } else if (dtype == 'dbinom') {
                                p.plus <- stats::pbinom(q = numReads['+'] - !right.is.bigger['+'], size = round(distr[as.character(states[,'pstate.right']), 'size'] * i1.bp/binsize), prob = distr[as.character(states[,'pstate.right']), 'prob'], lower.tail = right.is.bigger['+'])
                            }

                            p <- c(p.minus, p.plus)
                            p[is.na(p)] <- 1
                            p <- p['-'] * p['+']
                            ps[as.character(i1)] <- p
                        }
                        if (i1 >= 0) {
                            revps <- rev(ps)
                            if (revps[1] < 1-confint) {
                                i1 <- as.numeric(names(revps)[which(revps >= 1-confint)[1] - 1]) # first one from the end that is below threshold
                                if (is.na(i1)) {
                                    i1 <- as.numeric(names(ps)[1])
                                }
                            } else {
                                i1 <- as.numeric(names(revps)[1])
                            }
                            new.start <- end(cfrags)[ind-i1]
                            if (new.start < start(cbreaks.mod)[ibreak]) {
                                start(cbreaks.mod)[ibreak] <- new.start
                            }
                        }
                    }
                }
                
                ## Right side
                ind <- which(end(cfrags) > end(cbreaks)[ibreak])
                if (ibreak < length(cbreaks)) {
                    max.bp <- end(cbreaks[ibreak+1]) - end(cbreaks[ibreak])
                } else {
                    max.bp <- Inf
                }
                if (length(ind) > 0) {
                    ind <- ind[1]
                    p <- 1
                    ps <- numeric()
                    ps['0'] <- p
                    numReads <- c('-'=0, '+'=0)
                    if (!any(is.na(mus))) {
                        i1 <- -1
                        ## Helpers for speed improvement
                        end.cfrags <- end(cfrags)
                        end.cbreaks <- end(cbreaks)
                        strand.cfrags <- as.character(strand(cfrags))
                        while (p > 1-confint | i1 <= min.i1) {
                            i1 <- i1+1
                            if (ind+i1 > length(cfrags)) {
                                i1 = i1 - 1
                                break
                            }
                            i1.bp <- end.cfrags[ind+i1] - end.cbreaks[ibreak]
                            if (i1.bp > max.bp) {
                                i1 = i1 - 1
                                break
                            }
                            strandofread <- strand.cfrags[ind+i1]
                            numReads[strandofread] <- numReads[strandofread] + 1
                            
                            # Minus strand
                            dtype <- distr[as.character(states[,'mstate.left']), 'type']
                            if (dtype == 'dnbinom') {
                                p.minus <- stats::pnbinom(q = numReads['-'] - !left.is.bigger['-'], size = distr[as.character(states[,'mstate.left']), 'size'] * i1.bp/binsize, prob = distr[as.character(states[,'mstate.left']), 'prob'], lower.tail = left.is.bigger['-'])
                            } else if (dtype == 'dgeom') {
                                p.minus <- stats::pgeom(q = numReads['-'] - !left.is.bigger['-'], prob = dgeom.prob(distr[as.character(states[,'mstate.left']), 'mu'] * i1.bp/binsize), lower.tail = left.is.bigger['-'])
                            } else if (dtype == 'delta') {
                                p.minus <- c('-'=as.numeric(numReads['-'] == 0))
                            } else if (dtype == 'dbinom') {
                                p.minus <- stats::pbinom(q = numReads['-'] - !left.is.bigger['-'], size = round(distr[as.character(states[,'mstate.left']), 'size'] * i1.bp/binsize), prob = distr[as.character(states[,'mstate.left']), 'prob'], lower.tail = left.is.bigger['-'])
                            }
                            # Plus strand
                            dtype <- distr[as.character(states[,'pstate.left']), 'type']
                            if (dtype == 'dnbinom') {
                                p.plus <- stats::pnbinom(q = numReads['+'] - !left.is.bigger['+'], size = distr[as.character(states[,'pstate.left']), 'size'] * i1.bp/binsize, prob = distr[as.character(states[,'pstate.left']), 'prob'], lower.tail = left.is.bigger['+'])
                            } else if (dtype == 'dgeom') {
                                p.plus <- stats::pgeom(q = numReads['+'] - !left.is.bigger['+'], prob = dgeom.prob(distr[as.character(states[,'pstate.left']), 'mu'] * i1.bp/binsize), lower.tail = left.is.bigger['+'])
                            } else if (dtype == 'delta') {
                                p.plus <- c('+'=as.numeric(numReads['+'] == 0))
                            } else if (dtype == 'dbinom') {
                                p.plus <- stats::pbinom(q = numReads['+'] - !left.is.bigger['+'], size = round(distr[as.character(states[,'pstate.left']), 'size'] * i1.bp/binsize), prob = distr[as.character(states[,'pstate.left']), 'prob'], lower.tail = left.is.bigger['+'])
                            }
                            
                            p <- c(p.minus, p.plus)
                            p[is.na(p)] <- 1
                            p <- p['-'] * p['+']
                            ps[as.character(i1)] <- p
                        }
                        if (i1 >= 0) {
                            revps <- rev(ps)
                            if (revps[1] < 1-confint) {
                                i1 <- as.numeric(names(revps)[which(revps >= 1-confint)[1] - 1]) # first one from the end that is below threshold
                                if (is.na(i1)) {
                                    i1 <- as.numeric(names(ps)[1])
                                }
                            } else {
                                i1 <- as.numeric(names(revps)[1])
                            }
                            new.end <- start(cfrags)[ind+i1]
                            if (new.end > end(cbreaks.mod)[ibreak]) {
                                end(cbreaks.mod)[ibreak] <- new.end
                            }
                        }
                    }
                }
            }
        }
        breaks.conf[[chrom]] <- cbreaks.mod
    }
    breaks.conf <- unlist(breaks.conf, use.names = FALSE)
    stopTimedMessage(ptm)
    return(breaks.conf)
}
    

#' Refine breakpoints
#' 
#' Refine breakpoints with confidence intervals from an initial estimate (from \code{\link{getBreakpoints}}).
#' 
#' Breakpoints are refined by shifting the breakpoint within its initial confidence interval read by read and maximizing the probability of observing the left-right read distribution. 
#' 
#' @param model An \code{\link{aneuBiHMM}} object or a file that contains such an object.
#' @param breakpoints A \code{\link{GRanges}} object with breakpoints and confidence intervals, as returned by function \code{\link{getBreakpoints}}.
#' @param fragments A \code{\link{GRanges}} object with read fragments or a file that contains such an object.
#' @param confint Desired confidence interval for breakpoints.
#' @return A \code{\link{GRanges}} with breakpoint coordinates and confidence interals.
#' @export
#' 
#' @examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the data into bin size 1Mp
#'readfragments <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'), reads.return=TRUE)
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'## Fit the Hidden Markov Model
#'model <- findCNVs.strandseq(binned[[1]], eps=0.1, max.time=60)
#'## Add confidence intervals
#'breaks <- getBreakpoints(model, readfragments)
#'## Refine breakpoints
#'rfbreaks <- refineBreakpoints(model, breakpoints=breaks, fragments=readfragments)
#' 
refineBreakpoints <- function(model, breakpoints = model$breakpoints, fragments, confint=0.99) {
  
    ptm <- startTimedMessage("Refining breakpoints ...")
    model <- loadFromFiles(model, check.class = c("aneuHMM", "aneuBiHMM"))[[1]]
    fragments <- loadFromFiles(fragments, check.class = 'GRanges')[[1]]
    binsize <- mean(width(model$bincounts[[1]]))
    breaks <- breakpoints
    if (is.null(breaks)) {
        stop("No breakpoints found.")
    }
    if (is.null(breaks$start.conf)) {
        stop("No confidence intervals found.")
    }
    
    ## Sort fragments
    seqlevels(fragments) <- seqlevels(model$segments)
    fragments <- sort(fragments, ignore.strand=TRUE)
    fragments$p <- 0
    
    ## Distributions
    distr <- model$distributions
    if (class(distr) == 'list') {
        distr <- distr[[1]]
    }
    
    ## Do chromosomes one by one
    breaks.refined <- GenomicRanges::GRangesList()
    seqlevels(breaks.refined) <- seqlevels(breaks)
    for (chrom in unique(seqnames(breaks))) {
        # message("chrom = ", chrom)
        cbreaks <- breaks[seqnames(breaks) == chrom]
        cbreaks.mod <- cbreaks
        cfrags <- fragments[seqnames(fragments) == chrom]
        if (length(cbreaks) > 0) {
            for (ibreak in 1:length(cbreaks)) {
                # message(ibreak, "/", length(cbreaks))
              
                ## Go through confidence interval left to right
                states <- as.data.frame(mcols(cbreaks[ibreak]))
                ind <- which( (start(cfrags) >= cbreaks[ibreak]$start.conf) & (end(cfrags) <= cbreaks[ibreak]$end.conf) )
                frags.ind <- cfrags[ind]
                if (length(ind) > 0) {
                    numReads <- array(NA, dim=c(2,2,length(ind)+1), dimnames=list(strand=c('-','+'), direction=c('left','right'), ind=0:length(ind)))
                    ps <- array(NA, dim=c(length(ind)+1), dimnames=list(ind=0:length(ind)))
                    break.coords <- array(NA, dim=c(length(ind)+1), dimnames=list(ind=0:length(ind)))
                    ## Helpers for speed improvement
                    end.frags.ind <- end(frags.ind)
                    start.frags.ind <- start(frags.ind)
                    strand.frags.ind <- strand(frags.ind)
                    length.frags.ind <- length(frags.ind)
                    length.ind <- length(ind)
                    start.conf.cbreaks.ibreak <- cbreaks[ibreak]$start.conf
                    end.conf.cbreaks.ibreak <- cbreaks[ibreak]$end.conf
                    for (i1 in 0:length.ind) {
                        i1c <- as.character(i1)
                        if (i1 > 0 & i1 < length.ind) {
                            break.coord <- end.frags.ind[i1] + (start.frags.ind[i1+1] - end.frags.ind[i1]) / 2
                            numReads[, 'left', i1c] <- table(strand.frags.ind[1:i1])[c('-','+')]
                            numReads[, 'right', i1c] <- table(strand.frags.ind[(i1+1):length.frags.ind])[c('-','+')]
                        } else if (i1 == length.ind) {
                            break.coord <- end.conf.cbreaks.ibreak
                            numReads[, 'left', i1c] <- table(strand.frags.ind[1:i1])[c('-','+')]
                            numReads[, 'right', i1c] <- 0
                        } else {
                            break.coord <- start.conf.cbreaks.ibreak
                            numReads[, 'left', i1c] <- 0
                            numReads[, 'right', i1c] <- table(strand.frags.ind)[c('-','+')]
                        }
                        i1.bp <- break.coord - start.conf.cbreaks.ibreak

                        ptable <- numReads[,,1]
                        for (strand in dimnames(numReads)[[1]]) {
                            for (direction in dimnames(numReads)[[2]]) {
                                select <- paste0(c('-'='mstate.', '+'='pstate.')[strand], direction)
                                dtype <- distr[as.character(states[,select]), 'type']
                                if (dtype == 'dnbinom') {
                                    p <- dnbinom(x = numReads[strand, direction, i1c], size = distr[as.character(states[,select]), 'size'] * i1.bp/binsize, prob = distr[as.character(states[,select]), 'prob'])
                                } else if (dtype == 'dgeom') {
                                    p <- dgeom(x = numReads[strand, direction, i1c], prob = dgeom.prob(distr[as.character(states[,select]), 'mu'] * i1.bp/binsize))
                                } else if (dtype == 'delta') {
                                    p <- as.numeric(numReads[strand, direction, i1c] == 0)
                                } else if (dtype == 'dbinom') {
                                    p <- dnbinom(x = numReads[strand, direction, i1c], size = round(distr[as.character(states[,select]), 'size'] * i1.bp/binsize), prob = distr[as.character(states[,select]), 'prob'])
                                }
                                ptable[strand, direction] <- p
                            }
                        }
                        ps[i1c] <- prod(ptable)
                        break.coords[i1c] <- break.coord
                    }
                    new.break <- break.coords[names(ps)[which.max(ps)]]
                    if (new.break < start(cbreaks.mod)[ibreak]) {
                        start(cbreaks.mod)[ibreak] <- new.break
                        end(cbreaks.mod)[ibreak] <- new.break
                    } else {
                        end(cbreaks.mod)[ibreak] <- new.break
                        start(cbreaks.mod)[ibreak] <- new.break
                    }
                }
            }
        }
        breaks.refined[[chrom]] <- cbreaks.mod
    }
    breaks.refined <- unlist(breaks.refined, use.names = FALSE)
    stopTimedMessage(ptm)
    
    breaks.conf <- confidenceIntervals(breaks = breaks.refined, fragments = fragments, distr = distr, confint = confint, binsize = binsize)
    
    breaks.refined$start.conf <- start(breaks.conf)
    breaks.refined$end.conf <- end(breaks.conf)
    return(breaks.refined)
  
}



