

#' @import ggplot2
#' @import reshape2
#' @import grid
#' @import ggdendro
NULL

## Colors for plotting
#' aneufinder color scheme
#'
#' Get the color scheme that is used in the \pkg{\link{aneufinder}} plots.
#' @export
stateColors <- function() {
	state.colors <- c("mapped"="gray68","zero-inflation"="gray90", "nullsomy"="gray90","monosomy"="darkorchid2","disomy"="springgreen2","trisomy"="red3","tetrasomy"="gold2","multisomy"="deepskyblue","total"="black")
	return(state.colors)
}
#' aneufinder strand color scheme
#'
#' Get the color strand scheme that is used in the \pkg{\link{aneufinder}} plots.
#' @export
strandColors <- function() {
	return(c('+'="103,139,139", '-'="243,165,97"))
}

# =================================================================
# Define plotting methods for the generic
# =================================================================
#' Plotting function for saved \pkg{\link{aneufinder}} objects
#'
#' Convenience function that loads and plots a \pkg{\link{aneufinder}} object in one step.
#'
#' @param x A filename that contains either \code{\link{binned.data}} or a \code{\link{aneuHMM}}.
#' @param ... Additional arguments.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot character
#' @export
plot.character <- function(x, ...) {
	x <- get(load(x))
	plot(x, ...)
}

#' Plotting function for binned read counts
#'
#' Make plots for binned read counts from \code{\link{binned.data}}.
#'
#' @param x A \code{\link{GRanges}} object with binned read counts.
#' @param type Type of the plot, one of \code{c('karyogram', 'histogram', 'arrayCGH')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with read counts.}
#'   \item{\code{histogram}}{A histogram of read counts.}
#'   \item{\code{arrayCGH}}{An arrayCGH-like chromosome overview with read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRanges
#' @export
plot.GRanges <- function(x, type='karyogram', ...) {
	if (type == 'karyogram' | type==1) {
		plotKaryogram(x, ...)
	} else if (type == 'histogram' | type==2) {
		plotBinnedDataHistogram(x, ...)
	} else if (type == 'arrayCGH' | type==3) {
		plotArray(x, ...)
	}
}

#' Plotting function for \code{\link{aneuHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuHMM}} objects.
#'
#' @param x An \code{\link{aneuHMM}} object.
#' @param type Type of the plot, one of \code{c('karyogram', 'histogram', 'arrayCGH')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{An arrayCGH-like chromosome overview with CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuHMM
#' @export
plot.aneuHMM <- function(x, type='karyogram', ...) {
	if (type == 'karyogram' | type==1) {
		plotKaryogram(x, ...)
	} else if (type == 'histogram' | type==2) {
		plotUnivariateHistogram(x, ...)
	} else if (type == 'arrayCGH' | type==3) {
		plotArray(x, ...)
	}
}

#' Plotting function for \code{\link{aneuBiHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuBiHMM}} objects.
#'
#' @param x An \code{\link{aneuBiHMM}} object.
#' @param type Type of the plot, one of \code{c('karyogram', 'histogram', 'arrayCGH')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{An arrayCGH-like chromosome overview with CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuBiHMM
#' @export
plot.aneuBiHMM <- function(x, type='karyogram', ...) {
	if (type == 'karyogram' | type==1) {
		args <- names(list(...))
		if ('both.strands' %in% args) {
			plotKaryogram(x, ...)
		} else {
			plotKaryogram(x, both.strands=TRUE, ...)
		}
	} else if (type == 'histogram' | type==2) {
		plotBivariateHistograms(x, ...)
	} else if (type == 'arrayCGH' | type==3) {
		args <- names(list(...))
		if ('both.strands' %in% args) {
			plotArray(x, ...)
		} else {
			plotArray(x, both.strands=TRUE, ...)
		}
	}
}

# ============================================================
# Helper functions
# ============================================================
get_rightxlim <- function(counts) {
	rightxlim1 <- median(counts[counts>0])*7
	tab <- table(counts)
	tab <- tab[names(tab)!='0']
	breaks <- as.numeric(names(tab))
	rightxlim2 <- breaks[tab<=5 & breaks>median(counts)*2][1]
	rightxlim <- min(rightxlim1,rightxlim2, na.rm=TRUE)
	if (length(rightxlim)==0 | is.na(rightxlim) | is.infinite(rightxlim)) {
		rightxlim <- 1
	}
	return(rightxlim)
}

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges}} object.
#' @return The input \code{\link{GRanges}} with two additional metadata columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
	cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- seqlevels(gr)
	gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
	gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
	return(gr)
}

# ============================================================
# Plot a read histogram
# ============================================================
#' Plot a histogram of binned read counts
#'
#' Plot a histogram of binned read counts from \code{\link{binned.data}}
#'
#' @param binned.data A \code{\link{GRanges}} object containing binned read counts in meta-column 'counts'.
#' @param strand One of c('+','-','*'). Plot the histogram only for the specified strand.
#' @param chromosome,start,end Plot the histogram only for the specified chromosome, start and end position.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotBinnedDataHistogram <- function(binned.data, strand='*', chromosome=NULL, start=NULL, end=NULL) {

	# Select the rows to plot
	selectmask <- rep(TRUE,length(binned.data))
	numchrom <- length(table(seqnames(binned.data)))
	if (!is.null(chromosome)) {
		if (! chromosome %in% levels(seqnames(binned.data))) {
			stop(chromosome," can't be found in the binned data.")
		}
		selectchrom <- seqnames(binned.data) == chromosome
		selectmask <- selectmask & selectchrom
		numchrom <- 1
	}
	if (numchrom == 1) {
		if (!is.null(start)) {
			selectstart <- start(ranges(binned.data)) >= start
			selectmask <- selectmask & selectstart
		}
		if (!is.null(end)) {
			selectend <- end(ranges(binned.data)) <= end
			selectmask <- selectmask & selectend
		}
	}
	if (strand=='+') {
		select <- 'pcounts'
	} else if (strand=='-') {
		select <- 'mcounts'
	} else if (strand=='*') {
		select <- 'counts'
	}
	if (length(which(selectmask)) != length(binned.data$counts)) {
		counts <- mcols(binned.data)[,select][as.logical(selectmask)]
	} else {
		counts <- mcols(binned.data)[,select]
	}

	# Find the x limits
	breaks <- max(counts)
	if (max(counts)==0) { breaks <- 1 }
	rightxlim <- get_rightxlim(counts)

	# Plot the histogram
	ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")
	return(ggplt)

}

# =================================================================
# Plot a read histogram with univariate fits for a bivariate HMM
# =================================================================
plotBivariateHistograms <- function(bihmm) {

	binned.data <- bihmm$bins
	## Stack the two strands
	binned.data.minus <- binned.data
	strand(binned.data.minus) <- '-'
	binned.data.minus$counts <- binned.data.minus$mcounts
	binned.data.minus$counts.gc <- binned.data.minus$mcounts.gc
	binned.data.plus <- binned.data
	strand(binned.data.plus) <- '+'
	binned.data.plus$counts <- binned.data.plus$pcounts
	binned.data.plus$counts.gc <- binned.data.plus$pcounts.gc
	binned.data.stacked <- c(binned.data.minus, binned.data.plus)
	mask.attributes <- c(grep('complexity', names(attributes(binned.data)), value=T), 'spikyness', 'shannon.entropy')
	attributes(binned.data.stacked)[mask.attributes] <- attributes(binned.data)[mask.attributes]

	## Make fake uni.hmm and plot
	strand <- 'minus'
	uni.hmm <- list()
	uni.hmm$ID <- bihmm$ID
	uni.hmm$bins <- binned.data.stacked
	uni.hmm$bins$state <- uni.hmm$bins$pstate
	uni.hmm$bins$pstate <- NULL
	uni.hmm$bins$mstate <- NULL
	uni.hmm$weights <- bihmm$univariateParams$weights
	uni.hmm$distributions <- bihmm$distributions[[strand]]
	class(uni.hmm) <- class.univariate.hmm
	ggplts <- plotUnivariateHistogram(uni.hmm, strand='*')

# 	## Make fake uni.hmm and plot
# 	ggplts <- list()
# 	for (strand in c('plus','minus')) {
# 		uni.hmm <- list()
# 		uni.hmm$ID <- bihmm$ID
# 		uni.hmm$bins <- bihmm$bins
# 		uni.hmm$bins$state <- uni.hmm$bins$pstate
# 		uni.hmm$bins$pstate <- NULL
# 		uni.hmm$bins$mstate <- NULL
# 		uni.hmm$weights <- bihmm$weights.univariate[[strand]]
# 		uni.hmm$distributions <- bihmm$distributions[[strand]]
# 		class(uni.hmm) <- class.univariate.hmm
# 		ggplts[[strand]] <- plotUnivariateHistogram(uni.hmm, strand=strand)
# 	}
	
	return(ggplts)

}

# ============================================================
# Plot a read histogram with univariate fits
# ============================================================
#' Plot a histogram of binned read counts with fitted mixture distribution
#'
#' Plot a histogram of binned read counts from with fitted mixture distributions from a \code{\link{aneuHMM}} object.
#'
#' @param model A \code{\link{aneuHMM}} object.
#' @param state Plot the histogram only for the specified CNV-state.
#' @param strand One of c('+','-','*'). Plot the histogram only for the specified strand.
#' @param chromosome,start,end Plot the histogram only for the specified chromosome, start and end position.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotUnivariateHistogram <- function(model, state=NULL, strand='*', chromosome=NULL, start=NULL, end=NULL) {

	# Select the rows to plot
	selectmask <- rep(TRUE,length(model$bins))
	numchrom <- length(table(seqnames(model$bins)))
	if (!is.null(chromosome)) {
		if (! chromosome %in% levels(seqnames(model$bins))) {
			stop(chromosome," can't be found in the model coordinates.")
		}
		selectchrom <- seqnames(model$bins) == chromosome
		selectmask <- selectmask & selectchrom
		numchrom <- 1
	}
	if (numchrom == 1) {
		if (!is.null(start)) {
			selectstart <- start(ranges(model$bins)) >= start
			selectmask <- selectmask & selectstart
		}
		if (!is.null(end)) {
			selectend <- end(ranges(model$bins)) <= end
			selectmask <- selectmask & selectend
		}
	}
	if (!is.null(state)) {
		selectmask <- selectmask & model$bins$state==state
	}
	if (strand=='+' | strand=='plus') {
		select <- 'pcounts'
	} else if (strand=='-' | strand=='minus') {
		select <- 'mcounts'
	} else if (strand=='*' | strand=='both') {
		select <- 'counts'
	}
	if (length(which(selectmask)) != length(model$bins$counts)) {
		counts <- mcols(model$bins)[,select][as.logical(selectmask)]
	} else {
		counts <- mcols(model$bins)[,select]
	}
	states <- model$bins$state[as.logical(selectmask)]
	if (length(which(selectmask)) != length(model$bins)) {
		if (!is.null(state)) {
			weights <- rep(NA, length(levels(model$bins$state)))
			names(weights) <- levels(model$bins$state)
			for (istate in names(weights)) {
				weights[istate] <- length(which(states==levels(model$bins$state)[istate==names(weights)]))
			}
			weights <- weights / length(states)
		} else {
			weights <- model$weights
		}
	} else {
		if (!is.null(model$weights)) {
			weights <- model$weights
		}
	}

	# Find the x limits
	breaks <- max(counts)
	if (max(counts)==0) { breaks <- 1 }
	rightxlim <- get_rightxlim(counts)

	# Plot the histogram
	ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count") + ggtitle(model$ID)
	if (is.null(model$weights)) {
		return(ggplt)
	}

	### Add fits to the histogram
	c.state.labels <- as.character(levels(model$bins$state))
	numstates <- length(weights)
	x <- 0:max(counts)
	distributions <- data.frame(x)

	for (istate in 1:nrow(model$distributions)) {
		if (model$distributions[istate,'type']=='delta') {
			# zero-inflation
			distributions[[length(distributions)+1]] <- c(weights[istate],rep(0,length(x)-1))
		} else if (model$distributions[istate,'type']=='dgeom') {
			# geometric
			distributions[[length(distributions)+1]] <- weights[istate] * dgeom(x, model$distributions[istate,'prob'])
		} else if (model$distributions[istate,'type']=='dnbinom') {
			# negative binomials
			distributions[[length(distributions)+1]] <- weights[istate] * dnbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])
		} else if (model$distributions[istate,'type']=='dpois') {
			# poissons
			distributions[[length(distributions)+1]] <- weights[istate] * dpois(x, model$distributions[istate,'lambda'])
		} else if (model$distributions[istate,'type']=='dbinom') {
			# binomials
			s <- model$distributions[istate,'size']
			p <- model$distributions[istate,'prob']
# 			distributions[[length(distributions)+1]] <- weights[istate] * dbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])	# only defined for integer 'size'
			distributions[[length(distributions)+1]] <- weights[istate] * choose(s,x) * p^x * (1-p)^(s-x)
		}
	}
	distributions <- as.data.frame(distributions)
	names(distributions) <- c("x",c.state.labels)
	# Total
	distributions$total <- apply(distributions[-1], 1, sum)

	# Reshape the data.frame for plotting with ggplot
	distributions <- reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="state", times=c(c.state.labels,"total"))
	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions)
	} else {
		ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions[distributions$state==state,])
	}
	
	# Make legend and colors correct
	lmeans <- round(model$distributions[,'mu'], 2)
	lvars <- round(model$distributions[,'variance'], 2)
	legend <- paste(c.state.labels, ", mean=", lmeans, ", var=", lvars, sep='')
	legend <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=', round(var(counts),2)))
	ggplt <- ggplt + scale_color_manual(breaks=c(c.state.labels, 'total'), values=stateColors()[c(c.state.labels,'total')], labels=legend)
	ggplt <- ggplt + theme(legend.position=c(1,1), legend.justification=c(1,1))

	return(ggplt)

}


# ============================================================
# Plot karyogram-like chromosome overview
# ============================================================
#' Karyogram-like chromosome overview
#'
#' Plot a karyogram-like chromosome overview with read counts and CNV-state from a \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#'
#' @param model A \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#' @param file A PDF file where the plot will be saved.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotKaryogram <- function(model, both.strands=FALSE, file=NULL) {

	if (class(model)=='GRanges') {
		binned.data <- model
		model <- list()
		model$ID <- ''
		model$bins <- binned.data
		model$qualityInfo <- list(shannon.entropy=qc.entropy(binned.data$counts), spikyness=qc.spikyness(binned.data$counts), complexity=attr(binned.data, 'complexity.preseqR'), bhattacharyya=NA)
		plot.karyogram(model, both.strands=both.strands, file=file)
	} else if (class(model)==class.univariate.hmm) {
		plot.karyogram(model, both.strands=both.strands, file=file)
	} else if (class(model)==class.bivariate.hmm) {
		plot.karyogram(model, both.strands=both.strands, percentages=FALSE, file=file)
	}

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.karyogram <- function(model, both.strands=FALSE, percentages=TRUE, file=NULL) {
	
	## Convert to GRanges
	gr <- model$bins
	grl <- split(gr, seqnames(gr))

	## Get some variables
	num.chroms <- length(levels(seqnames(gr)))
	maxseqlength <- max(seqlengths(gr))
	tab <- table(gr$counts)
	tab <- tab[names(tab)!='0']
	custom.xlim1 <- as.numeric(names(tab)[which.max(tab)]) # maximum value of read distribution
	custom.xlim2 <- as.integer(mean(gr$counts, trim=0.05)) # mean number of counts
	custom.xlim <- max(custom.xlim1, custom.xlim2, na.rm=TRUE) * 4
	if (both.strands) {
		custom.xlim <- custom.xlim / 2
	}
	if (length(custom.xlim)==0 | is.na(custom.xlim)) {
		custom.xlim <- 1
	}

	## Setup page
	fs_title <- 20
	fs_x <- 13
	nrows <- 2	# rows for plotting chromosomes
	nrows.text <- 2	# two additional rows for displaying ID and qualityInfo
	nrows.total <- nrows + nrows.text
	ncols <- ceiling(num.chroms/nrows)
	if (!is.null(file)) {
		pdf(file=file, width=ncols*1.4, height=nrows*4.6)
	}
	grid.newpage()
	layout <- matrix(1:((nrows.total)*ncols), ncol=ncols, nrow=nrows.total, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,1,21,21))))
	# Main title
	grid.text(model$ID, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_title))
	# Quality info
	if (is.null(model$qualityInfo$complexity)) { model$qualityInfo$complexity <- NA }
	if (is.null(model$qualityInfo$spikyness)) { model$qualityInfo$spikyness <- NA }
	if (is.null(model$qualityInfo$shannon.entropy)) { model$qualityInfo$shannon.entropy <- NA }
	if (is.null(model$qualityInfo$bhattacharyya)) { model$qualityInfo$bhattacharyya <- NA }
	quality.string <- paste0('complexity = ',round(model$qualityInfo$complexity),',  spikyness = ',round(model$qualityInfo$spikyness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2),',  bhattacharyya = ',round(model$qualityInfo$bhattacharyya,2))
	grid.text(quality.string, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_x))

	## Get SCE coordinates
	if (both.strands) {
		scecoords <- model$sce
		# Set to midpoint
		start(scecoords) <- (start(scecoords)+end(scecoords))/2
		end(scecoords) <- start(scecoords)
	}

	## Go through chromosomes and plot
	for (i1 in 1:num.chroms) {
		# Get the i,j matrix positions of the regions that contain this subplot
		matchidx <- as.data.frame(which(layout == i1+nrows.text*ncols, arr.ind = TRUE))
		if (!is.null(grl[[i1]]$state) & percentages) {
			# Percentage of chromosome in state
			tstate <- table(mcols(grl[[i1]])$state)
			pstate.all <- tstate / sum(tstate)
			pstate <- round(pstate.all*100)[-1]	# without 'nullsomy / unmapped' state
			pstring <- apply(pstate, 1, function(x) { paste0(": ", x, "%") })
			pstring <- paste0(names(pstring), pstring)
			pstring <- paste(pstring[which.max(pstate)], collapse="\n")
			pstring2 <- round(pstate.all*100)[1]	# only 'nullsomy / unmapped'
			pstring2 <- paste0(names(pstring2), ": ", pstring2, "%")
		} else {
			pstring <- ''
			pstring2 <- ''
		}

		# Plot the read counts
		dfplot <- as.data.frame(grl[[i1]])
		# Add offset so that end coordinates match (to plot p-arm on top)
			offset <- (maxseqlength-seqlengths(gr)[i1])
			dfplot$start <- dfplot$start + offset
			dfplot$end <- dfplot$end + offset
		# Set values too big for plotting to limit
			dfplot$counts[dfplot$counts>=custom.xlim] <- custom.xlim
			dfplot.points <- dfplot[dfplot$counts>=custom.xlim,]
			dfplot.points$counts <- rep(custom.xlim, nrow(dfplot.points))

			if (both.strands) {
				dfplot$mcounts <- - dfplot$mcounts	# negative minus counts
				dfplot$pcounts[dfplot$pcounts>=custom.xlim] <- custom.xlim
				dfplot$mcounts[dfplot$mcounts<=-custom.xlim] <- -custom.xlim
				dfplot.points.plus <- dfplot[dfplot$pcounts>=custom.xlim,]
				dfplot.points.plus$counts <- rep(custom.xlim, nrow(dfplot.points.plus))
				dfplot.points.minus <- dfplot[dfplot$mcounts<=-custom.xlim,]
				dfplot.points.minus$counts <- rep(-custom.xlim, nrow(dfplot.points.minus))
			}

		empty_theme <- theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
			axis.title.x=element_text(size=fs_x),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
		ggplt <- ggplot(dfplot, aes_string(x='start', y='counts'))	# data
		if (!is.null(grl[[i1]]$state)) {
			if (both.strands) {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts', col='pstate'), size=0.2)	# read count
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts', col='mstate'), size=0.2)	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts', col='pstate'), size=5, shape=21)	# outliers
				ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts', col='mstate'), size=5, shape=21)	# outliers
			} else {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts', col='state'), size=0.2)	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts', col='state'), size=2, shape=21)	# outliers
			}
			ggplt <- ggplt + scale_color_manual(values=stateColors(), drop=F)	# do not drop levels if not present
		} else {
			if (both.strands) {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts'), size=0.2, col='gray20')	# read count
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts'), size=0.2, col='gray20')	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')	# outliers
				ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')	# outliers
			} else {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts'), size=0.2, col='gray20')	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts'), size=2, shape=21, col='gray20')	# outliers
			}
		}
		if (both.strands) {
			ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=-maxseqlength, xmax=-offset, col='white', fill='gray20')	# chromosome backbone as simple rectangle, the minus sign in x-axis is a dirty hack because of a bug in ggplt when reversing axis
		} else {
			ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=-maxseqlength, xmax=-offset, col='white', fill='gray20')	# chromosome backbone as simple rectangle
		}
		if (both.strands) {
			dfsce <- as.data.frame(scecoords[seqnames(scecoords)==names(grl)[i1]])
			dfsce$start <- dfsce$start + offset
			dfsce$end <- dfsce$end + offset
			if (nrow(dfsce)>0) {
				ggplt <- ggplt + geom_segment(data=dfsce, aes(x=start, xend=start), y=-custom.xlim, yend=-0.5*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
			}
		}
		ggplt <- ggplt + empty_theme	# no axes whatsoever
		ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1], "\n", pstring, "\n", pstring2))	# chromosome names
		if (both.strands) {
			ggplt <- ggplt + xlim(maxseqlength,0) + ylim(-custom.xlim,custom.xlim)	# set x- and y-limits
		} else {
			ggplt <- ggplt + xlim(maxseqlength,0) + ylim(-0.6*custom.xlim,custom.xlim)	# set x- and y-limits
		}
		ggplt <- ggplt + coord_flip() # flip coordinates
		suppressWarnings(
		print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		)
		
	}
	if (!is.null(file)) {
		d <- dev.off()
	}
}


# =================================================================
# Plot a heatmap of chromosome state for multiple samples
# =================================================================
#' Plot aneuploidy state
#'
#' Plot a heatmap of aneuploidy state for multiple samples. Samples can be clustered and the output can be returned as data.frame.
#'
#' @param hmm.list A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmm.list}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity in their CNV-state.
#' @param as.data.frame If \code{TRUE}, instead of a plot, a data.frame with the aneuploidy state for each sample will be returned.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or a data.frame, depending on option \code{as.data.frame}.
#' @author Aaron Taudt
#' @export
heatmapAneuploidies <- function(hmm.list, ylabels=NULL, cluster=TRUE, as.data.frame=FALSE) {

	## Check user input
	if (!is.null(ylabels)) {
		if (length(ylabels) != length(hmm.list)) {
			stop("length(ylabels) must equal length(hmm.list)")
		}
	}

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)
	levels.state <- unique(unlist(lapply(hmm.list, function(hmm) { levels(hmm$bins$state) })))
	
	## Assign new IDs
	if (!is.null(ylabels)) {
		for (i1 in 1:length(hmm.list)) {
			hmm.list[[i1]]$ID <- ylabels[i1]
		}
	}

	## Transform to GRanges in reduced representation
	grlred <- GRangesList()
	for (hmm in hmm.list) {
		if (!is.null(hmm$segments)) {
			grlred[[hmm$ID]] <- hmm$segments
		}
	}
	
	## Find the most frequent state (mfs) for each chromosome and sample
	message("finding most frequent state for each sample and chromosome ...", appendLF=F); ptm <- proc.time()
	grl.per.chrom <- lapply(grlred, function(x) { split(x, seqnames(x)) })
	mfs.samples <- list()
	for (i1 in 1:length(grlred)) {
		mfs.samples[[names(grlred)[i1]]] <- lapply(grl.per.chrom[[i1]], function(x) {
      if (length(x)>0) {
        tab <- aggregate(width(x), by=list(state=x$state), FUN="sum")
        tab$state[which.max(tab$x)]
      } else {
        "nullsomy"
      }
      })
		attr(mfs.samples[[names(grlred)[i1]]], "varname") <- 'chromosome'
	}
	attr(mfs.samples, "varname") <- 'sample'
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Transform to data.frame
	# Long format
	df <- reshape2::melt(mfs.samples, value.name='state')
	df$state <- factor(df$state, levels=levels.state)
	df$sample <- factor(df$sample, levels=unique(df$sample))
	df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
	# Wide format
	df.wide <- reshape2::dcast(df, sample ~ chromosome, value.var='state', factorsAsStrings=F)
	# Correct strings to factors
	for (col in 2:ncol(df.wide)) {
		df.wide[,col] <- factor(df.wide[,col], levels=levels.state)
	}

	## Cluster the samples by chromosome state
	if (cluster) {
		# Cluster
		hc <- hclust(dist(data.matrix(df.wide[-1])))
		# Reorder samples in mfs list
		mfs.samples.clustered <- mfs.samples[hc$order]
		attr(mfs.samples.clustered, "varname") <- 'sample'
		df <- reshape2::melt(mfs.samples.clustered, value.name='state')
		df$state <- factor(df$state, levels=levels.state)
		df$sample <- factor(df$sample, levels=unique(df$sample))
		df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
	}

	## Plot to heatmap
	if (as.data.frame) {
		df.table <- df.wide
		for (i1 in 2:ncol(df.table)) {
			df.table[,i1] <- as.numeric(df.table[,i1])-1
		}
		return(df.table)
	} else {
		ggplt <- ggplot(df) + geom_tile(aes_string(x='chromosome', y='sample', fill='state'), col='black') + theme_bw() + scale_fill_manual(values=stateColors()[levels(df$state)])
		return(ggplt)
	}
}


# =================================================================
# Plot a clustered heatmap of state calls
# =================================================================
#' Genome wide heatmap of CNV-state
#'
#' Plot a genome wide heatmap of copy number variation state. This heatmap is best plotted to file, because in most cases it will be too big for cleanly plotting it to screen.
#'
#' @param hmm.list A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmm.list}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param plot.SCE Logical indicating whether SCE events should be plotted.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @export
heatmapGenomewide <- function(hmm.list, ylabels=NULL, file=NULL, cluster=TRUE, plot.SCE=TRUE) {

	## Check user input
	if (!is.null(ylabels)) {
		if (length(ylabels) != length(hmm.list)) {
			stop("length(ylabels) must equal length(hmm.list)")
		}
	}

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Assign new IDs
	if (!is.null(ylabels)) {
		for (i1 in 1:length(hmm.list)) {
			hmm.list[[i1]]$ID <- ylabels[i1]
		}
	}

	## Get segments and SCE coordinates
	temp <- getSegments(hmm.list, cluster=cluster)
	grlred <- temp$segments
	hc <- temp$clustering
	if (cluster) {
		hmm.list <- hmm.list[hc$order]
	}
	if (plot.SCE) {
		sce <- lapply(hmm.list,'[[','sce')
		sce <- sce[!unlist(lapply(sce, is.null))]
		sce <- sce[lapply(sce, length)!=0]
		if (length(sce)==0) {
			plot.SCE <- FALSE
		}
	}

	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
	message("transforming coordinates ...", appendLF=F); ptm <- proc.time()
	grlred <- endoapply(grlred, transCoord)
	if (plot.SCE) {
		sce <- endoapply(sce, transCoord)
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Data.frame for plotting
	message("Making the plot ...", appendLF=F); ptm <- proc.time()
	# Data
	df <- list()
	for (i1 in 1:length(grlred)) {
		df[[length(df)+1]] <- data.frame(start=grlred[[i1]]$start.genome, end=grlred[[i1]]$end.genome, seqnames=seqnames(grlred[[i1]]), sample=names(grlred)[i1], state=grlred[[i1]]$state)
	}
	df <- do.call(rbind, df)
	if (plot.SCE) {
		df.sce <- list()
		for (i1 in 1:length(sce)) {
			df.sce[[length(df.sce)+1]] <- data.frame(start=sce[[i1]]$start.genome, end=sce[[i1]]$end.genome, seqnames=seqnames(sce[[i1]]), sample=names(grlred)[i1])
		}
		df.sce <- do.call(rbind, df.sce)
	}
	# Chromosome lines
	cum.seqlengths <- cumsum(as.numeric(seqlengths(grlred[[1]])))
	names(cum.seqlengths) <- seqlevels(grlred[[1]])
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- seqlevels(grlred[[1]])
	label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grlred[[1]]) )
	df.chroms <- data.frame(y=c(0,cum.seqlengths))

	### Plot ###

	## Prepare the plot
	ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=stateColors()) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20))
	ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
	if (plot.SCE) {
		ggplt <- ggplt + geom_point(data=df.sce, mapping=aes_string(x='sample', y='start'), size=2) + ylab('')
	}
	## Prepare the dendrogram
	if (!is.null(hc)) {
		dhc <- as.dendrogram(hc)
		ddata <- dendro_data(dhc, type = "rectangle")
		ggdndr <- ggplot(segment(ddata)) + geom_segment(aes_string(x='x', y='y', xend='xend', yend='yend')) + coord_flip() + scale_y_reverse(expand=c(0,0)) + theme_dendro()
	}

	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	printToGrid <- function() {
# 		grob.dndr <- ggplotGrob(ggdndr)
# 		grob.plt <- ggplotGrob(ggplt)
# 		grob.dndr <- gtable_add_cols(grob.dndr, unit(width.cm,'cm'))
# 		grob.dndr <- gtable_add_grob(grob.dndr, h, t=1, l=ncol(d), b=3, r=ncol(d))
# 		grid.newpage()
# 		grid.draw(grob.dndr)
		# Setup page
		grid.newpage()
		layout <- matrix(1:2, ncol=2, nrow=1, byrow=T)
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), widths=c(10,90))))
		# Print to page
		print(ggplt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
		if (!is.null(hc)) {
			print(ggdndr, vp=viewport(layout.pos.row=1, layout.pos.col=1))
		}
	}

	## Plot to file
	width.cm <- sum(as.numeric(seqlengths(hmm.list[[1]]$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
	if (!is.null(file)) {
		message("plotting to file ",file," ...", appendLF=F); ptm <- proc.time()
		height.cm <- length(hmm.list) * 0.5
		width.dendro.cm <- 20
		pdf(file, width=(width.cm+width.dendro.cm)/2.54, height=height.cm/2.54)
		printToGrid()
# 		pdf(file, width=width.cm/2.54, height=height.cm/2.54)
# 		print(ggplt)
		d <- dev.off()
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	} else {
# 		printToGrid()
		return(ggplt)
	}

}


# =================================================================
# Plot arrayCGH-like chromosome overview
# =================================================================
#' ArrayCGH-like chromosome overview
#'
#' Plot an arrayCGH-like chromosome overview with read counts and CNV-state from a \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#'
#' @param model A \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#' @param file A PDF file where the plot will be saved.
#' @param plot.SCE Logical indicating whether SCE events should be plotted.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotArray <- function(model, both.strands=FALSE, plot.SCE=TRUE, file=NULL) {

	if (class(model)=='GRanges') {
		binned.data <- model
		model <- list()
		model$ID <- ''
		model$bins <- binned.data
		model$qualityInfo <- list(shannon.entropy=qc.entropy(binned.data$counts), spikyness=qc.spikyness(binned.data$counts), complexity=attr(binned.data, 'complexity.preseqR'))
		plot.array(model, both.strands=both.strands, plot.SCE=FALSE, file=file)
	} else if (class(model)==class.univariate.hmm) {
		plot.array(model, both.strands=FALSE, plot.SCE=FALSE, file=file)
	} else if (class(model)==class.bivariate.hmm) {
		plot.array(model, both.strands=both.strands, plot.SCE=plot.SCE, file=file)
	}

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.array <- function(model, both.strands=FALSE, plot.SCE=TRUE, file=NULL) {
	
	## Convert to GRanges
	bins <- model$bins
	## Get SCE coordinates
	if (plot.SCE) {
		scecoords <- model$sce
		# Set to midpoint
		start(scecoords) <- (start(scecoords)+end(scecoords))/2
		end(scecoords) <- start(scecoords)
	}

	## Get some variables
	num.chroms <- length(levels(seqnames(bins)))
	maxseqlength <- max(seqlengths(bins))
	tab <- table(bins$counts)
	tab <- tab[names(tab)!='0']
	custom.xlim1 <- as.numeric(names(tab)[which.max(tab)]) # maximum value of read distribution
	custom.xlim2 <- as.integer(mean(bins$counts, trim=0.05)) # mean number of counts
	custom.xlim <- max(custom.xlim1, custom.xlim2, na.rm=TRUE) * 2.7
	if (both.strands) {
		custom.xlim <- custom.xlim / 1
	}
	if (length(custom.xlim)==0 | is.na(custom.xlim)) {
		custom.xlim <- 1
	}

	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
	message("transforming coordinates ...", appendLF=F); ptm <- proc.time()
	bins <- transCoord(bins)
	if (plot.SCE) {
		scecoords <- transCoord(scecoords)
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Setup page
	fs_title <- 20
	fs_x <- 13
	nrows <- 1	# rows for plotting chromosomes
	nrows.text <- 2	# two additional rows for displaying ID and qualityInfo
	nrows.total <- nrows + nrows.text
	ncols <- 1
	if (!is.null(file)) {
		pdf(file=file, width=20, height=5)
	}
	grid.newpage()
	layout <- matrix(1:((nrows.total)*ncols), ncol=ncols, nrow=nrows.total, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,1,21))))
	# Main title
	grid.text(model$ID, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_title))
	# Quality info
	if (is.null(model$qualityInfo$complexity)) { model$qualityInfo$complexity <- NA }
	if (is.null(model$qualityInfo$spikyness)) { model$qualityInfo$spikyness <- NA }
	if (is.null(model$qualityInfo$shannon.entropy)) { model$qualityInfo$shannon.entropy <- NA }
	if (is.null(model$qualityInfo$bhattacharyya)) { model$qualityInfo$bhattacharyya <- NA }
	quality.string <- paste0('complexity = ',round(model$qualityInfo$complexity),',  spikyness = ',round(model$qualityInfo$spikyness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2),',  bhattacharyya = ',round(model$qualityInfo$bhattacharyya,2))
	grid.text(quality.string, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_x))

	# Get the i,j matrix positions of the regions that contain this subplot
	matchidx <- as.data.frame(which(layout == 1+nrows.text*ncols, arr.ind = TRUE))

	# Plot the read counts
	dfplot <- as.data.frame(bins)
	# Set values too big for plotting to limit
	if (both.strands) {
		dfplot$mcounts <- - dfplot$mcounts	# negative minus counts
	}
	# Mean counts for CNV-state
	if (!is.null(model$segments$state)) {
		dfplot.seg <- as.data.frame(transCoord(model$segments))
		if (class(model)==class.univariate.hmm) {
			dfplot.seg$counts.CNV <- model$distributions[as.character(dfplot.seg$state),'mu']
		} else if (class(model)==class.bivariate.hmm) {
			dfplot.seg$counts.CNV <- model$distributions$plus[as.character(dfplot.seg$state),'mu']
			dfplot.seg$pcounts.CNV <- model$distributions$plus[as.character(dfplot.seg$pstate),'mu']
			dfplot.seg$mcounts.CNV <- -model$distributions$minus[as.character(dfplot.seg$mstate),'mu']
		}
	}

	empty_theme <- theme(axis.line=element_blank(),
		axis.text.x=element_text(size=fs_x),
		axis.ticks=element_blank(),
		axis.title.x=element_blank(),
		panel.background=element_blank(),
		panel.border=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		plot.background=element_blank())
	ggplt <- ggplot(dfplot, aes_string(x='start.genome', y='counts'))	# data
	if (both.strands) {
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='pcounts'), position=position_jitter(width=0))	# read count
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='mcounts'), position=position_jitter(width=0))	# read count
	} else {
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='counts'), position=position_jitter(width=0))	# read count
	}
	if (!is.null(bins$state)) {
		if (both.strands) {
			ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='pcounts.CNV',xend='end.genome',yend='pcounts.CNV', col='pstate'), size=2)
			ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='mcounts.CNV',xend='end.genome',yend='mcounts.CNV', col='mstate'), size=2)
		} else {
			ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='counts.CNV',xend='end.genome',yend='counts.CNV', col='state'), size=2)
		}
	}
	# Chromosome lines
	cum.seqlengths <- cumsum(as.numeric(seqlengths(bins)))
	names(cum.seqlengths) <- seqlevels(bins)
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- seqlevels(bins)
	label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(bins) )
	df.chroms <- data.frame(x=c(0,cum.seqlengths))
	ggplt <- ggplt + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2)
	
	ggplt <- ggplt + scale_color_manual(values=stateColors(), drop=F)	# do not drop levels if not present
	if (plot.SCE) {
		dfsce <- as.data.frame(scecoords)
		if (nrow(dfsce)>0) {
			ggplt <- ggplt + geom_segment(data=dfsce, aes_string(x='start.genome', xend='start.genome'), y=-1.5*custom.xlim, yend=-1.3*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
		}
	}
	ggplt <- ggplt + empty_theme	# no axes whatsoever
	if (both.strands) {
		ggplt <- ggplt + ylim(-1.5*custom.xlim,custom.xlim)	# set x- and y-limits
	} else {
		ggplt <- ggplt + ylim(0,custom.xlim)	# set x- and y-limits
	}
	# Get midpoints of each chromosome for xticks
	ggplt <- ggplt + scale_x_continuous(breaks=seqlengths(model$bins)/2+cum.seqlengths.0[as.character(seqlevels(model$bins))], labels=seqlevels(model$bins))
	ggplt <- ggplt + ylab('read count')
	suppressWarnings(
	print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
	)
		
	if (!is.null(file)) {
		d <- dev.off()
	}
}

