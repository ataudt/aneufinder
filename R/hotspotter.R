#' Find hotspots of genomic events
#'
#' Find hotspots of genomic events by using kernel \link{density} estimation.
#'
#' The hotspotter uses \code{\link[stats]{density}} to perform a KDE. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events (bootstrapping).
#' 
#' @param breakpoints A list with \code{\link{GRanges-class}} object containing the coordinates of the genomic events.
#' @param bw Bandwidth used for kernel density estimation (see \code{\link[stats]{density}}).
#' @param pval P-value cutoff for hotspots.
#' @param spacing.bp Spacing of datapoints for KDE in basepairs.
#' @return A list of \code{\link{GRanges-class}} objects containing 1) coordinates of hotspots and 2) p-values within the hotspot.
#' @importFrom stats ecdf p.adjust runif
#' @importFrom S4Vectors endoapply
#' @author Aaron Taudt
hotspotter <- function(breakpoints, bw, pval=5e-2, spacing.bp=5000) {

  ptm <- startTimedMessage("Finding breakpoint hotspots ...")
  set.seed(0) # fix seed for random permutations of bootstrapping
  
	## Coerce into one GRanges
	names(breakpoints) <- NULL
	gr <- do.call(c,breakpoints)
	gr <- sort(gr)
	
	## Iterate over chromosomes and calculate p-values
	pranges.list <- GRangesList()
	pvalues.list <- GRangesList()
	for (chrom in seqlevels(gr)) {
		grc <- gr[seqnames(gr)==chrom]
		if (length(grc)>1) {
			midpoints <- (start(grc)+end(grc))/2
  		n <- seqlengths(gr)[chrom] / spacing.bp
			kde <- stats::density(midpoints, bw=bw, kernel='gaussian', n=n)
			# Random distribution of genomic events
			kde.densities <- numeric()
			for (i1 in 1:100) {
  			midpoints.r <- round(stats::runif(length(midpoints),1,seqlengths(gr)[chrom]))
  			kde.r <- stats::density(midpoints.r, bw=bw, kernel='gaussian', n=n)
  			kde.densities <- c(kde.densities, kde.r$y)
			}
			# Use ecdf to calculate p-values 
			p <- 1-stats::ecdf(kde.densities)(kde$y)
			p <- p.adjust(p = p, method = 'BY')
			pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
			# Make GRanges
			pvalues$end <- pvalues$start
			pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
			pvalues$kde <- kde$y
			pvalues <- as(pvalues,'GRanges')
			seqlevels(pvalues) <- seqlevels(gr)
			suppressWarnings(
				seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
			)
			# # Resize from pointsize to spacing
			# suppressWarnings(
			# 	pvalues <- resize(pvalues, width=spacing.bp, fix='center')
			# )
			pvalues <- trim(pvalues)
			## Find regions where p-value is below specification
			mask <- pvalues$pvalue <= pval
			rle.pvals <- rle(mask)
			rle.pvals$values <- cumsum(rle.pvals$values+1)
			pvalues$group <- inverse.rle(rle.pvals)
			pvalues.list[[chrom]] <- pvalues[mask]
			if (length(which(mask))>0) {
				pvalues.split <- split(pvalues[mask],pvalues$group[mask])
				pranges.split <- GRangesList()
				for (i1 in 1:length(pvalues.split)) {
				  ipvalues <- pvalues.split[[i1]]
  				max.density <- max(ipvalues$kde)
  				min.pvalue <- min(ipvalues$pvalue)
  				ipvalues.max <- ipvalues[ipvalues$kde==max.density]
  				pranges <- GRanges(seqnames=chrom, ranges=IRanges(start=start(ipvalues)[1], end=end(ipvalues)[length(ipvalues)]), strand='*')
  				seqlevels(pranges) <- seqlevels(gr)
  				pranges$start.max <- start(ipvalues.max)[1]
  				pranges$end.max <- end(ipvalues.max)[length(ipvalues.max)]
  				pranges$pval <- min.pvalue
  				pranges.split[[i1]] <- pranges
				}
				pranges <- unlist(pranges.split, use.names = FALSE)
      	pranges <- resize(pranges, fix='center', width=width(pranges) + spacing.bp)
				pranges$num.events <- countOverlaps(pranges, grc)
				pranges.list[[chrom]] <- pranges
			}
		}
	}
	pranges <- unlist(pranges.list, use.names=FALSE)
	pranges <- trim(pranges)
	pvalues <- unlist(pvalues.list, use.names=FALSE)
	pvalues <- trim(pvalues)
	pvalues$group <- NULL
	
	stopTimedMessage(ptm)
	return(list(hotspots = pranges, densities = pvalues))

}


#' Find hotspots of genomic events
#'
#' Find hotspots of genomic events by using kernel density estimation.
#'
#' The hotspotter uses a gaussian kernel with variable bandwidth to perform a KDE. The bandwidth depends on the confidence intervals of the breakpoints. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events (bootstrapping). 
#' 
#' @param breakpoints A list with \code{\link{GRanges-class}} object containing the coordinates of the genomic events and their confidence intervals.
#' @param confint Confidence interval that was used for breakpoint estimation.
#' @param pval P-value cutoff for hotspots.
#' @param spacing.bp Spacing of datapoints for KDE in basepairs.
#' @return A list of \code{\link{GRanges-class}} objects containing 1) coordinates of hotspots and 2) p-values within the hotspot.
#' @importFrom stats ecdf p.adjust runif dnorm
#' @importFrom S4Vectors endoapply
#' @importFrom utils head tail
#' @author Aaron Taudt
hotspotter.variable <- function(breakpoints, confint, pval=5e-2, spacing.bp=5000) {

    ptm <- startTimedMessage("Finding breakpoint hotspots ...")
    set.seed(0) # fix seed for random permutations of bootstrapping
    ## Function for rolling sum
    rsum.cumsum <- function(x, n) {
        x <- utils::tail(cumsum(x) - cumsum(c(rep(0, n), utils::head(x, -n))), -n + 1)
        if (n %% 2 == 0) {
            x <- c(rep(0, n/2), x, rep(0, n/2-1))
        } else {
            x <- c(rep(0, n/2), x, rep(0, n/2))
        }
        return(x)
    }
    
  	## Coerce into one GRanges
  	names(breakpoints) <- NULL
  	gr <- do.call(c,breakpoints)
  	gr <- sort(gr)
  	midpoints <- (start(gr)+end(gr))/2
  	ranges(gr) <- IRanges(start=gr$start.conf, end=gr$end.conf)
  	mcols(gr) <- NULL
  	gr$midpoints <- midpoints
  	
  	## Iterate over chromosomes and calculate p-values
  	pranges.list <- GRangesList()
  	pvalues.list <- GRangesList()
  	for (chrom in seqlevels(gr)) {
    		grc <- gr[seqnames(gr)==chrom]
    		
      	# Make spaced datapoints for KDE
        s <- seq(from=1, to=seqlengths(gr)[chrom], by=spacing.bp)
        kdecoord <- GRanges(seqnames=chrom, ranges=IRanges(start=s, end=s))
        
    		if (length(grc)>1) {
    		  
            ## Go through breakpoints and add densities
            # bandwidths <- sqrt(-(width(grc)/2)^2/2/log(1-confint))
            left.bandwidths <- sqrt(-(grc$midpoints - start(grc))^2/2/log(1-confint))
            right.bandwidths <- sqrt(-(end(grc) - grc$midpoints)^2/2/log(1-confint))
            kde <- rep(0, length(s))
            for (i1 in 1:length(grc)) {
                # kde <- kde + stats::dnorm(s, mean=grc$midpoints[i1], sd = bandwidths[i1])
                midpoint <- grc$midpoints[i1]
                mask <- s <= midpoint
                kde[mask] <- kde[mask] + stats::dnorm(s[mask], mean=midpoint, sd = left.bandwidths[i1])
                kde[!mask] <- kde[!mask] + stats::dnorm(s[!mask], mean=midpoint, sd = right.bandwidths[i1])
            }
            # KDE integral over median bandwidth
            # median.bw <- median(bandwidths)
            median.bw <- median(left.bandwidths + right.bandwidths)
            points2sum <- round(median.bw / spacing.bp)
            kde <- rsum.cumsum(kde, points2sum)
            
            ## Repeat for randomly reshuffled breakpoints (bootstrapping)
      			kde.densities <- numeric()
      			for (ibootstrap in 1:10) {
          			midpoints.r <- round(stats::runif(n=length(grc), min=1, max=seqlengths(gr)[chrom]))
                # Go through breakpoints and add densities
                kde.bootstrap <- rep(0, length(s))
                for (i1 in 1:length(grc)) {
                    # kde.bootstrap <- kde.bootstrap + stats::dnorm(s, mean=midpoints.r[i1], sd = bandwidths[i1])
                    midpoint <- midpoints.r[i1]
                    mask <- s <= midpoint
                    kde.bootstrap[mask] <- kde.bootstrap[mask] + stats::dnorm(s[mask], mean=midpoint, sd = left.bandwidths[i1])
                    kde.bootstrap[!mask] <- kde.bootstrap[!mask] + stats::dnorm(s[!mask], mean=midpoint, sd = right.bandwidths[i1])
                }
          			# kde.densities <- c(kde.densities, kde.bootstrap)
                kde.densities <- c(kde.densities, rsum.cumsum(kde.bootstrap, points2sum))
      			}
    	
      			# Use ecdf to calculate p-values 
      			p <- 1-stats::ecdf(kde.densities)(kde)
      			p <- p.adjust(p = p, method = 'BY')
      			pvalues <- data.frame(chromosome=chrom,start=s,pvalue=p)
      			# Make GRanges
      			pvalues$end <- pvalues$start
      			pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
      			pvalues$kde <- kde
      			pvalues <- as(pvalues,'GRanges')
      			seqlevels(pvalues) <- seqlevels(gr)
      			suppressWarnings(
        				seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
      			)
      			pvalues <- trim(pvalues)
      			## Find regions where p-value is below specification
      			mask <- pvalues$pvalue <= pval
      			rle.pvals <- rle(mask)
      			rle.pvals$values <- cumsum(rle.pvals$values+1)
      			pvalues$group <- inverse.rle(rle.pvals)
      			pvalues.list[[chrom]] <- pvalues[mask]
      			if (length(which(mask))>0) {
        				pvalues.split <- split(pvalues[mask],pvalues$group[mask])
        				pranges.split <- GRangesList()
        				for (i1 in 1:length(pvalues.split)) {
          				  ipvalues <- pvalues.split[[i1]]
            				max.density <- max(ipvalues$kde)
            				min.pvalue <- min(ipvalues$pvalue)
            				ipvalues.max <- ipvalues[ipvalues$kde==max.density]
            				pranges <- GRanges(seqnames=chrom, ranges=IRanges(start=start(ipvalues)[1], end=end(ipvalues)[length(ipvalues)]), strand='*')
            				seqlevels(pranges) <- seqlevels(gr)
            				pranges$start.max <- start(ipvalues.max)[1]
            				pranges$end.max <- end(ipvalues.max)[length(ipvalues.max)]
            				pranges$pval <- min.pvalue
            				pranges.split[[i1]] <- pranges
        				}
        				pranges <- unlist(pranges.split, use.names = FALSE)
              	pranges <- resize(pranges, fix='center', width=width(pranges) + spacing.bp)
        				pranges$num.events <- countOverlaps(pranges, grc)
        				pranges.list[[chrom]] <- pranges
      			}
    		}
  	}
  	pranges <- unlist(pranges.list, use.names=FALSE)
  	pranges <- trim(pranges)
  	pvalues <- unlist(pvalues.list, use.names=FALSE)
  	pvalues <- trim(pvalues)
  	pvalues$group <- NULL
	
  	stopTimedMessage(ptm)
  	return(list(hotspots = pranges, densities = pvalues))
}


#' Find breakpoint hotspots
#' 
#' Find breakpoint hotspots with kernel density estimation (KDE).
#' 
#' \code{findHotspots} uses \code{\link[stats]{density}} to perform a KDE. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events. Due to this random sampling, the result can vary for each function call, most likely for hotspots whose p-value is close to the specified \code{pval}.
#' 
#' @param models A list of \code{\link{GRanges-class}} or \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @inheritParams hotspotter
#' @param filename Will write hotspot coordinates and densities to the specified file. Endings "_breakpoint-hotspots.bed.gz" and "_breakpoint-densities.wig.gz" will be appended to \code{filename}.
#' @return A list of \code{\link{GRanges-class}} objects containing 1) coordinates of hotspots and 2) p-values within the hotspot.
#' @export
#' 
findHotspots <- function(models, bw, pval=5e-2, spacing.bp=5000, filename=NULL) {
  
    models <- loadFromFiles(models, check.class=c("aneuHMM", "aneuBiHMM"))
    
    ## Extract breakpoints
    ptm <- startTimedMessage("Extracting breakpoints ...")
    breakpoints <- list()
    total.read.count <- numeric()
    for (i1 in 1:length(models)) {
        if (is.character(models[[i1]])) {
            file <- models[[i1]]
            hmm <- suppressMessages( loadFromFiles(file)[[1]] )
        } else {
            file <- names(models)[i1]
            hmm <- models[[i1]]
        }
        breakpoints[[basename(file)]] <- hmm$breakpoints
        total.read.count[basename(file)] <- hmm$qualityInfo$total.read.count
    }
    if (is.null(bw)) {
        bw <- sum(as.numeric(seqlengths(hmm$bins))) / mean(total.read.count)
    }
    stopTimedMessage(ptm)
    ptm <- startTimedMessage("Estimating hotspots ...")
    hslist <- hotspotter(breakpoints, bw=bw, pval=pval, spacing.bp=spacing.bp)
    stopTimedMessage(ptm)
    
    ## Writing to file
    if (!is.null(filename)) {
        ## Breakpoint hotspots
        savename <- paste0(filename, '_breakpoint-hotspots')
        if (!file.exists(paste0(savename,'.bed.gz'))) {
            hotspots <- hslist$hotspots
            if (!is.null(hotspots)) {
                exportGRanges(hotspots, filename=savename, trackname=basename(savename), score=hotspots$num.events, thickStart = hotspots$start.max, thickEnd = hotspots$end.max, priority=41)
            }
        }
        ## Hotspot densities
        savename <- paste0(filename, '_breakpoint-hotspot-densities')
        if (!file.exists(paste0(savename,'.wig.gz'))) {
            densities <- hslist$densities
            if (!is.null(hotspots)) {
                exportGRanges(densities, filename=savename, trackname=basename(savename), as.wiggle = TRUE, wiggle.val = densities$kde, priority=40)
            }
        }
      
    }
    return(hslist)
}
