#' Find hotspots of genomic events
#'
#' Find hotspots of genomic events by using kernel \link{density} estimation.
#'
#' @param gr.list A list with \code{\link{GRanges}} object containing the coordinates of the genomic events.
#' @param bw Bandwidth used for kernel density estimation (see \code{\link[stats]{density}}).
#' @param pval P-value cutoff for hotspots.
#' @return A \code{\link{GRanges}} object containing coordinates of hotspots with p-values.
#' @author Aaron Taudt
#' @export
hotspotter <- function(gr.list, bw, pval=0.05) {

	## Coerce into one GRanges
	names(gr.list) <- NULL
	gr <- do.call(c,gr.list)
	gr <- sort(gr)
	
	## Iterate over chromosomes and calculate p-values
	pvalues <- list()
	for (chrom in seqlevels(gr)) {
		grc <- gr[seqnames(gr)==chrom]
		midpoints <- (start(grc)+end(grc))/2
		kde <- stats::density(midpoints,bw=bw,kernel='gaussian')
		# Random distribution of genomic events
		midpoints.r <- round(runif(midpoints,1,seqlengths(gr)[chrom]))
		kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian')
		# Use ecdf to calculate p-values 
		p <- 1-ecdf(kde.r$y)(kde$y)
		pvalues[[chrom]] <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
	}
	pvalues <- do.call(rbind,pvalues)
	rownames(pvalues) <- NULL
	pvalues$end <- pvalues$start
	pvalues <- as(pvalues,'GRanges')
	suppressWarnings(
		seqlengths(pvalues) <- seqlengths(gr.list[[1]])
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
	pvalues.split <- split(pvalues[mask],pvalues$group[mask])
	pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
	pranges$group <- NULL
	pranges$num.events <- countOverlaps(pranges,gr)

	return(pranges)

}
