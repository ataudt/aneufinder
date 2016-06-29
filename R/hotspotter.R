#' Find hotspots of genomic events
#'
#' Find hotspots of genomic events by using kernel \link{density} estimation.
#'
#' The hotspotter uses \code{\link[stat]{density}} to perform a KDE. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events. Due to this random sampling, the result can vary for each function call, most likely for hotspots whose p-value is close to the specified \code{pval}.
#' 
#' @param gr.list A list with \code{\link{GRanges}} object containing the coordinates of the genomic events.
#' @param bw Bandwidth used for kernel density estimation (see \code{\link[stats]{density}}).
#' @param pval P-value cutoff for hotspots.
#' @return A \code{\link{GRanges}} object containing coordinates of hotspots with p-values.
#' @importFrom stats ecdf p.adjust runif
#' @author Aaron Taudt
hotspotter <- function(gr.list, bw, pval=1e-8) {

	## Coerce into one GRanges
	names(gr.list) <- NULL
	gr <- do.call(c,gr.list)
	gr <- sort(gr)
	
	## Iterate over chromosomes and calculate p-values
	pranges.list <- GRangesList()
	for (chrom in seqlevels(gr)) {
		grc <- gr[seqnames(gr)==chrom]
		if (length(grc)>1) {
			midpoints <- (start(grc)+end(grc))/2
			kde <- stats::density(midpoints,bw=bw,kernel='gaussian')
			# Random distribution of genomic events
			kde.densities <- numeric()
			for (i1 in 1:100) {
  			midpoints.r <- round(stats::runif(length(midpoints),1,seqlengths(gr)[chrom]))
  			kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian')
  			kde.densities <- c(kde.densities, kde.r$y)
			}
			# Use ecdf to calculate p-values 
			p <- 1-stats::ecdf(kde.densities)(kde$y)
			pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
			# Make GRanges
			pvalues$end <- pvalues$start
			pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
			pvalues <- as(pvalues,'GRanges')
			suppressWarnings(
				seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
			)
			# Resize from pointsize to bandwidth
			suppressWarnings(
				pvalues <- resize(pvalues, width=bw, fix='center')
			)
			pvalues <- trim(pvalues)
			## Find regions where p-value is below specification
			mask <- pvalues$pvalue <= pval
			rle.pvals <- rle(mask)
			rle.pvals$values <- cumsum(rle.pvals$values+1)
			pvalues$group <- inverse.rle(rle.pvals)
			if (length(which(mask))>0) {
				pvalues.split <- split(pvalues[mask],pvalues$group[mask])
				pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
				pranges$group <- NULL
				pranges$num.events <- countOverlaps(pranges,grc)
				pranges.list[[chrom]] <- pranges
			}
		}
	}
	pranges <- unlist(pranges.list, use.names=FALSE)
	names(pranges) <- NULL

	return(pranges)

}
