#' Find hotspots of genomic events
#'
#' Find hotspots of genomic events by using kernel \link{density} estimation.
#'
#' The hotspotter uses \code{\link[stats]{density}} to perform a KDE. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events. Due to this random sampling, the result can vary for each function call, most likely for hotspots whose p-value is close to the specified \code{pval}.
#' 
#' @param gr.list A list with \code{\link{GRanges}} object containing the coordinates of the genomic events.
#' @param bw Bandwidth used for kernel density estimation (see \code{\link[stats]{density}}).
#' @param pval P-value cutoff for hotspots.
#' @param spacing.bp Spacing of datapoints for KDE in basepairs.
#' @return A list of \code{\link{GRanges}} objects containing 1) coordinates of hotspots and 2) p-values within the hotspot.
#' @importFrom stats ecdf p.adjust runif
#' @author Aaron Taudt
hotspotter <- function(gr.list, bw, pval=5e-2, spacing.bp=5000) {

	## Coerce into one GRanges
	names(gr.list) <- NULL
	gr <- do.call(c,gr.list)
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
	pvalues <- unlist(pvalues.list, use.names=FALSE)
	pvalues$group <- NULL
	

	return(list(hotspots = pranges, densities = pvalues))

}
