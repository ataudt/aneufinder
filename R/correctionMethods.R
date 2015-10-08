

#' GC correction
#'
#' Correct a list of \code{\link{binned.data}} by GC content
#'
#' @param binned.data.list A \code{list()} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
#' @param same.GC.content If \code{TRUE} the GC content will only be calculated once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects describe the same genome at the same binsize.
#' @author Aaron Taudt
#' @importFrom Biostrings Views alphabetFrequency
#' @export
correctGC <- function(binned.data.list, GC.BSgenome, same.GC.content=FALSE) {

	binned.data.list <- loadGRangesFromFiles(binned.data.list)
	same.GC.calculated <- FALSE
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]

		## Check if seqlengths of data and GC.correction are consistent
		# Replace 1->chr1 if necessary
		chromlengths <- seqlengths(binned.data)
		chroms <- names(chromlengths)
		mask <- !grepl('^chr', chroms)
		chroms[mask] <- paste0('chr',chroms[mask])
		names(chromlengths) <- chroms
		# Compare
		compare <- chromlengths[chroms] == seqlengths(GC.BSgenome)[chroms]
		if (any(compare==FALSE, na.rm=T)) {
			warning(paste0(attr(binned.data,'ID'),": Chromosome lengths differ between binned data and 'GC.BSgenome'. GC correction skipped. Please use the correct genome for option 'GC.BSgenome'."))
			binned.data.list[[i1]] <- binned.data
			next
		}

		## Calculate GC content per bin
		if (same.GC.content & !same.GC.calculated | !same.GC.content) {
			message("Calculating GC content per bin ...", appendLF=F); ptm <- proc.time()
			GC.content <- list()
			for (chrom in seqlevels(binned.data)) {
				if (!grepl('^chr',chrom)) {
					chr <- paste0('chr',chrom)
				} else {
					chr <- chrom
				}
				if (chr %in% seqlevels(GC.BSgenome)) {
					view <- Biostrings::Views(GC.BSgenome[[chr]], ranges(binned.data)[seqnames(binned.data)==chrom])
					freq <- Biostrings::alphabetFrequency(view, as.prob = T, baseOnly=T)
					if (nrow(freq) > 1) {
						GC.content[[as.character(chrom)]] <- rowSums(freq[, c("G","C")])
					} else {
						GC.content[[as.character(chrom)]] <- sum(freq[, c("G","C")])
					}
				} else {
					GC.content[[as.character(chrom)]] <- rep(NA, length(binned.data[seqnames(binned.data)==chrom]))
					warning(paste0(attr(binned.data,'ID'),": No sequence information for chromosome ",chr," available."))
				}
			}
			GC.content <- unlist(GC.content)
			same.GC.calculated <- TRUE
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		binned.data$GC <- GC.content

		### GC correction ###
		message("GC correction ...", appendLF=F); ptm <- proc.time()
		counts <- binned.data$counts
		mcounts <- binned.data$mcounts
		pcounts <- binned.data$pcounts
		## Correction factors
		gc.categories <- seq(from=0, to=1, length=20)
		intervals.per.bin <- findInterval(binned.data$GC, gc.categories)
		intervals <- sort(unique(intervals.per.bin))
		mean.counts.global <- mean(binned.data$counts, trim=0.05)
		correction.factors <- NULL
		weights <- NULL
		for (interval in intervals) {
			mask <- intervals.per.bin==interval
			counts.with.same.GC <- binned.data$counts[mask]
			weights[as.character(gc.categories[interval])] <- length(counts.with.same.GC)
			mean.counts.with.same.GC <- mean(counts.with.same.GC, na.rm=T, trim=0.05)
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
		fit <- lm(y ~ poly(x, 2, raw=T), data=df, weights=weight)
		fitted.correction.factors <- predict(fit, data.frame(x=gc.categories[intervals]))
		names(fitted.correction.factors) <- intervals
		for (interval in intervals) {
			mask <- which(intervals.per.bin==interval)
			correction.factor <- fitted.correction.factors[as.character(interval)]
			counts[mask] <- counts[mask] * correction.factor
			mcounts[mask] <- mcounts[mask] * correction.factor
			pcounts[mask] <- pcounts[mask] * correction.factor
		}
		binned.data$counts <- as.integer(round(counts))
		binned.data$mcounts <- as.integer(round(mcounts))
		binned.data$pcounts <- as.integer(round(pcounts))
		# Produce fit to check
		ggplt <- ggplot(df) + geom_point(aes_string(x='x', y='y', size='weight')) + geom_line(aes_string(x='x', y='y'), data=data.frame(x=gc.categories[intervals], y=fitted.correction.factors)) + theme_bw() + ggtitle('GC correction') + xlab('GC content') + ylab('correction factor')
		attr(binned.data, 'GC.correction.ggplt') <- ggplt
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		### Quality measures ###
		## Spikyness
		attr(binned.data, 'spikyness') <- qc.spikyness(binned.data$counts)
		## Shannon entropy
		attr(binned.data, 'shannon.entropy') <- qc.entropy(binned.data$counts)

		binned.data.list[[i1]] <- binned.data
	}
	return(binned.data.list)
}


#' Mappability correction
#'
#' Correct a list of \code{\link{binned.data}} by mappability. WARNING: This approach works only for a set of super heterogeneous samples.
#'
#' @param binned.data.list A \code{list()} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @author Aaron Taudt
#' @export
correctMappability <- function(binned.data.list) {

	bins <- loadGRangesFromFiles(binned.data.list)
	if (length(bins)<=1) {
		stop("argument 'hmms' expects a list of length(hmms)>=2")
	}
	## Check if all models have the same binsize
	binsizes <- unlist(lapply(bins, function(bin) { width(bin)[1] }))
	if (any(binsizes!=binsizes[1])) {
		warning("No mappability correction done. All binsizes have to be the same.")
		return(bins)
	}
	## Variables
	num.bins <- length(bins[[1]])
	## Get all counts
	counts <- matrix(NA, nrow=length(bins[[1]]), ncol=length(bins))
	for (i1 in 1:length(bins)) {
		counts[,i1] <- bins[[i1]]$counts
	}
	## Normalize to the mean number of counts per bin
	ncounts <- sweep(counts, 2, colSums(counts), '/') * num.bins * mean(counts)
	mcounts <- apply(ncounts, 1, mean)
	## Correction factor
	cfactor <- mean(mcounts) / mcounts
	cfactor[is.infinite(cfactor)] <- 0
	for (i1 in 1:length(bins)) {
		bins[[i1]]$counts.map <- as.integer(round(bins[[i1]]$counts * cfactor))
		bins[[i1]]$mcounts.map <- as.integer(round(bins[[i1]]$mcounts * cfactor))
		bins[[i1]]$pcounts.map <- as.integer(round(bins[[i1]]$pcounts * cfactor))
	}
	return(bins)
}
