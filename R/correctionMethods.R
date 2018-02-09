#' GC correction
#'
#' Correct a list of \code{\link{binned.data}} by GC content.
#'
#' Two methods are available for GC correction:  Option \code{method='quadratic'} uses the method described in the Supplementary of \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit to adjust the read count.
#' 
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
#' @param same.binsize If \code{TRUE} the GC content will only be calculated once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects describe the same genome at the same binsize and stepsize.
#' @param method One of \code{c('quadratic', 'loess')}. Option \code{method='quadratic'} uses the method described in the Supplementary of \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit to adjust the read count.
#' @param return.plot Set to \code{TRUE} if plots should be returned for visual assessment of the GC correction.
#' @param bins A \code{\link{binned.data}} object with meta-data column 'GC'. If this is specified, \code{GC.BSgenome} is ignored. Beware, no format checking is done.
#' @return A \code{list()} with \code{\link{binned.data}} objects with adjusted read counts. Alternatively a \code{list()} with \code{\link[ggplot2]{ggplot}} objects if \code{return.plot=TRUE}.
#' @author Aaron Taudt
#' @importFrom Biostrings Views alphabetFrequency
#' @importFrom stats lm predict loess
#' @importFrom reshape2 melt
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
correctGC <- function(binned.data.list, GC.BSgenome, same.binsize=FALSE, method='loess', return.plot=FALSE, bins=NULL) {

  if (is.null(bins)) {
      ## Determine format of GC.BSgenome
      if (grepl('^chr', seqlevels(GC.BSgenome)[1])) {
        bsgenome.format <- 'UCSC'
      } else {
        bsgenome.format <- 'NCBI'
      }
  } else {
      ## Determine format of bins
      if (grepl('^chr', seqlevels(bins)[1])) {
        bsgenome.format <- 'UCSC'
      } else {
        bsgenome.format <- 'NCBI'
      }
  }
  
  ### Loop over all bin entries ###
	binned.data.list <- loadFromFiles(binned.data.list, check.class=c('GRanges','GRangesList'))
	same.binsize.calculated <- FALSE
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]

		if (is.null(bins)) {
  		## Check if seqlengths of data and GC.correction are consistent
  		# Replace 1->chr1 if necessary
  		chromlengths <- seqlengths(binned.data)
  		chroms <- names(chromlengths)
  		if (bsgenome.format == 'UCSC') {
    		mask <- !grepl('^chr', chroms)
    		chroms[mask] <- paste0('chr',chroms[mask])
  		} else if (bsgenome.format == 'NCBI') {
    		mask <- grepl('^chr', chroms)
    		chroms[mask] <- sub('^chr', '', chroms[mask])
  		}
  		names(chromlengths) <- chroms
  		# Compare
  		compare <- chromlengths[chroms] == seqlengths(GC.BSgenome)[chroms]
  		if (any(compare==FALSE, na.rm=TRUE)) {
  			warning(paste0(attr(binned.data,'ID'),": Incorrect 'GC.BSgenome' specified. seqlengths() differ. GC correction skipped. Please use the correct genome for option 'GC.BSgenome'."))
  			binned.data.list[[i1]] <- binned.data
  			next
  		}
		}

		if (class(binned.data) == 'GRanges') {
		  blist <- GRangesList('0'=binned.data)
		  attr(blist, 'qualityInfo') <- attr(binned.data, 'qualityInfo')
		} else if (is(binned.data, "GRangesList")) {
		  blist <- binned.data
		}
		## Calculate GC content per bin
		if (same.binsize & !same.binsize.calculated | !same.binsize) {
  			ptm <- startTimedMessage("Calculating GC content per bin and stepsize ...")
    		GC <- list()
  			for (i2 in 1:length(blist)) {
    			  iblist <- blist[[i2]]
      			GC.content <- list()
      			for (chrom in seqlevels(iblist)) {
        				if (!grepl('^chr', chrom) & bsgenome.format == 'UCSC') {
          					chr <- paste0('chr', chrom)
        				} else if (grepl('^chr', chrom) & bsgenome.format == 'NCBI') {
          				  chr <- sub('^chr', '', chrom)
        				} else {
          					chr <- chrom
        				}
        			  if (is.null(bins)) {
            				if (chr %in% seqlevels(GC.BSgenome)) {
            					view <- Biostrings::Views(GC.BSgenome[[chr]], ranges(iblist)[seqnames(iblist)==chrom])
            					freq <- Biostrings::alphabetFrequency(view, as.prob = TRUE, baseOnly=TRUE)
            					GC.content[[as.character(chrom)]] <- rowSums(freq[, c("G","C"), drop=FALSE])
            				} else {
            					GC.content[[as.character(chrom)]] <- rep(NA, length(iblist[seqnames(iblist)==chrom]))
            					warning(paste0(attr(iblist,'ID'),": No sequence information for chromosome ",chr," available."))
            				}
        			  } else {
        			      GC.content[[as.character(chrom)]] <- bins[[i2]][seqnames(bins[[i2]])==chrom]$GC
        			  }
      			}
      			GC.content <- unlist(GC.content)
      			GC[[i2]] <- GC.content
  			}
  			same.binsize.calculated <- TRUE
  			stopTimedMessage(ptm)
		}
		
		### GC correction ###
		ptm <- startTimedMessage("GC correction for sample ", i1, " ...")
		blist.gc <- GRangesList()
		plots <- list()
		for (i2 in 1:length(blist)) {
		  iblist <- blist[[i2]]
		  iblist$GC <- GC[[i2]]

  		counts <- iblist$counts
  		mcounts <- iblist$mcounts
  		pcounts <- iblist$pcounts
  		if (method == 'quadratic') {
    		## Correction factors
    		gc.categories <- seq(from=0, to=1, length=20)
    		intervals.per.bin <- findInterval(iblist$GC, gc.categories)
    		intervals <- sort(unique(intervals.per.bin))
    		mean.counts.global <- mean(iblist$counts, trim=0.05)
    		correction.factors <- NULL
    		weights <- NULL
    		for (interval in intervals) {
    			mask <- intervals.per.bin==interval
    			counts.with.same.GC <- iblist$counts[mask]
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
    		if (return.plot) {
      		# Produce fit to check
      		# ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y', size='weight')) + geom_line(aes_string(x='x', y='y'), col='red', data=data.frame(x=gc.categories[intervals], y=fitted.correction.factors)) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('correction factor')
    		  df <- data.frame(counts=iblist$counts, GC=iblist$GC)
    		  df$counts.GC <- counts
    		  dfplot <- df[,c('counts','GC','counts.GC')]
    		  dfplot <- reshape2::melt(dfplot, id.vars='GC', variable.name='method', value.name='counts')
    		  dfplot$method <- c(counts='counts', counts.GC='GC-corrected counts')[dfplot$method]
    		  ggplt <- ggplot(dfplot) + geom_point(aes_string(x='GC', y='counts'), alpha=0.1) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('counts') + facet_wrap(~ method)
    		  # Mean counts
    		  ggplt <- ggplt + geom_hline(yintercept = mean.counts.global, linetype=2)
    		  # Add quadratic fit
    		  dffit <- data.frame(GC=df$GC, fit=mean.counts.global / stats::predict(fit, data.frame(x=gc.categories[intervals.per.bin])), fit.GC=NA)
    		  dfplot <- reshape2::melt(dffit, id.vars='GC', variable.name='method', value.name='counts')
    		  dfplot$method <- c(fit='counts', fit.GC='GC-corrected counts')[dfplot$method]
    		  ggplt <- ggplt + geom_line(mapping = aes_string(x='GC', y='counts'), data=dfplot, col='red')
      		plots[[i2]] <- ggplt
      		ggplt + coord_cartesian(xlim=c(0,1))
    		}
  		
  		} else if (method == 'loess') {
    		mean.counts.global <- mean(iblist$counts, trim=0.05)
  		  df <- as.data.frame(mcols(iblist))
  		  fit <- stats::loess(counts ~ GC, data=df)
  		  correction.factor <- mean.counts.global / fit$fitted
  			counts <- counts * correction.factor
  			mcounts <- mcounts * correction.factor
  			pcounts <- pcounts * correction.factor
    		if (return.plot) {
      		# Produce fit to check
    # 		  df$counts.scaled <- df$counts / mean.counts.global
    # 		  df$correction.factor <- correction.factor
    #   		ggplt <- ggplot(df) + geom_point(aes_string(x='GC', y='counts.scaled'), alpha=0.1) + geom_line(aes_string(x='GC', y='correction.factor'), col='red') + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('correction factor')
    		  df$counts.GC <- counts
    		  dfplot <- df[,c('counts','GC','counts.GC')]
    		  dfplot <- reshape2::melt(dfplot, id.vars='GC', variable.name='method', value.name='counts')
    		  dfplot$method <- c(counts='counts', counts.GC='GC-corrected counts')[dfplot$method]
    		  ggplt <- ggplot(dfplot) + geom_point(aes_string(x='GC', y='counts'), alpha=0.1) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('counts') + facet_wrap(~ method)
    		  # Mean counts
    		  ggplt <- ggplt + geom_hline(yintercept = mean.counts.global, linetype=2)
    		  # Add loess fit
    		  dffit <- data.frame(GC=df$GC, fit=fit$fitted, fit.GC=NA)
    		  dfplot <- reshape2::melt(dffit, id.vars='GC', variable.name='method', value.name='counts')
    		  dfplot$method <- c(fit='counts', fit.GC='GC-corrected counts')[dfplot$method]
    		  ggplt <- ggplt + geom_line(mapping = aes_string(x='GC', y='counts'), data=dfplot, col='red')
      		plots[[i2]] <- ggplt
    		}
  		}
  		iblist$counts <- as.integer(round(counts))
  		iblist$mcounts <- as.integer(round(mcounts))
  		iblist$pcounts <- as.integer(round(pcounts))
  		iblist$counts[iblist$counts<0] <- 0
  		iblist$pcounts[iblist$pcounts<0] <- 0
  		iblist$mcounts[iblist$mcounts<0] <- 0
  
		  blist.gc[[i2]] <- iblist
		}
		names(blist.gc) <- names(blist)
		attr(blist.gc, 'qualityInfo') <- attr(blist, 'qualityInfo')
		stopTimedMessage(ptm)
		
		if (!return.plot) {
  		## Spikyness
  		attr(blist.gc, 'qualityInfo')$spikiness <- qc.spikiness(blist.gc[[1]]$counts)
  		## Shannon entropy
  		attr(blist.gc, 'qualityInfo')$entropy <- qc.entropy(blist.gc[[1]]$counts)
  		## ID
  		attr(blist.gc, 'ID') <- attr(binned.data, 'ID')
  		# Return
  		binned.data.list[[i1]] <- blist.gc
		} else {
		  binned.data.list[[i1]] <- plots
		}

	}
	return(binned.data.list)
}


