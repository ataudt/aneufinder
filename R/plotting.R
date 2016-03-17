

#' @import ggplot2
#' @importFrom cowplot plot_grid draw_label
#' @import reshape2
#' @import ggdendro
NULL

#' \pkg{AneuFinder} color scheme
#'
#' Get the color schemes that are used in the AneuFinder plots.
#'
#' @return A character vector with colors.
#' @name colors
#' @aliases stateColors strandColors
NULL

#' @describeIn colors Colors that are used for the states.
#' @param states A character vector with states whose color should be returned.
#' @export
#'@examples
#'## Make a nice pie chart with the AneuFinder state color scheme
#'statecolors <- stateColors()
#'pie(rep(1,length(statecolors)), labels=names(statecolors), col=statecolors)
#'
stateColors <- function(states=c('zero-inflation', paste0(0:10, '-somy'), 'total')) {
	state.colors <- c("zero-inflation"="gray90", "0-somy"="gray90","1-somy"="darkorchid2","2-somy"="springgreen2","3-somy"="red3","4-somy"="gold2","5-somy"="lightpink4","6-somy"="lightpink3","7-somy"="lightpink2","8-somy"="lightpink1","9-somy"="lightpink","10-somy"="deepskyblue","total"="black")
	states.with.color <- intersect(states, names(state.colors))
	cols <- rep('black', length(states))
	names(cols) <- states
	cols[states.with.color] <- state.colors[states.with.color]
	return(cols)
}

#' @describeIn colors Colors that are used to distinguish strands.
#' @param strands A character vector with strands whose color should be returned. Any combination of \code{c('+','-','*')}.
#' @export
#'@examples
#'## Make a nice pie chart with the AneuFinder strand color scheme
#'strandcolors <- strandColors()
#'pie(rep(1,length(strandcolors)), labels=names(strandcolors), col=strandcolors)
#'
strandColors <- function(strands=c('+','-')) {
	strand.colors <- c('+'="#678B8B", '-'="#F3A561", '*'="#000000")
	strands.with.color <- intersect(strands, names(strand.colors))
	cols <- rep('black', length(strands))
	names(cols) <- strands
	cols[strands.with.color] <- strand.colors[strands.with.color]
	return(cols)
}

# =================================================================
# Define plotting methods for the generic
# =================================================================
#' Plotting function for saved \pkg{\link{AneuFinder}} objects
#'
#' Convenience function that loads and plots a \pkg{\link{AneuFinder}} object in one step.
#'
#' @param x A filename that contains either \code{\link{binned.data}} or a \code{\link{aneuHMM}}.
#' @param ... Additional arguments.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot character
#' @importFrom graphics plot
#' @export
plot.character <- function(x, ...) {
	x <- get(load(x))
	graphics::plot(x, ...)
}

#' Plotting function for binned read counts
#'
#' Make plots for binned read counts from \code{\link{binned.data}}.
#'
#' @param x A \code{\link{GRanges}} object with binned read counts.
#' @param type Type of the plot, one of \code{c('profile', 'histogram', 'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with read counts.}
#'   \item{\code{histogram}}{A histogram of read counts.}
#'   \item{\code{profile}}{An profile with read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRanges
#' @export
plot.GRanges <- function(x, type='profile', ...) {
	if (type == 'karyogram' | type==3) {
		plotKaryogram(x, ...)
	} else if (type == 'histogram' | type==2) {
		plotBinnedDataHistogram(x, ...)
	} else if (type == 'profile' | type==1) {
		plotProfile(x, ...)
	}
}

#' Plotting function for \code{\link{aneuHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuHMM}} objects.
#'
#' @param x An \code{\link{aneuHMM}} object.
#' @param type Type of the plot, one of \code{c('profile', 'histogram', 'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{An profile with read counts and CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuHMM
#' @export
plot.aneuHMM <- function(x, type='profile', ...) {
	if (type == 'karyogram' | type==3) {
		plotKaryogram(x, ...)
	} else if (type == 'histogram' | type==2) {
		plotUnivariateHistogram(x, ...)
	} else if (type == 'profile' | type==1) {
		plotProfile(x, ...)
	}
}

#' Plotting function for \code{\link{aneuBiHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuBiHMM}} objects.
#'
#' @param x An \code{\link{aneuBiHMM}} object.
#' @param type Type of the plot, one of \code{c('profile', 'histogram', 'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{profile}}{An profile with read counts and CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuBiHMM
#' @export
plot.aneuBiHMM <- function(x, type='profile', ...) {
	if (type == 'karyogram' | type==3) {
		args <- names(list(...))
		if ('both.strands' %in% args) {
			plotKaryogram(x, ...)
		} else {
			plotKaryogram(x, both.strands=TRUE, ...)
		}
	} else if (type == 'histogram' | type==2) {
		plotBivariateHistograms(x, ...)
	} else if (type == 'profile' | type==1) {
		args <- names(list(...))
		if ('both.strands' %in% args) {
			plotProfile(x, ...)
		} else {
			plotProfile(x, both.strands=TRUE, ...)
		}
	}
}

# ============================================================
# Helper functions
# ============================================================
get_rightxlim <- function(counts) {
# 	rightxlim1 <- median(counts[counts>0])*7
# 	tab <- table(counts)
# 	tab <- tab[names(tab)!='0']
# 	breaks <- as.numeric(names(tab))
# 	rightxlim2 <- breaks[tab<=5 & breaks>median(counts)*2][1]
# 	rightxlim <- min(rightxlim1,rightxlim2, na.rm=TRUE)
	rightxlim <- stats::quantile(counts, 0.999)
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
	mask.attributes <- c(grep('complexity', names(attributes(binned.data)), value=TRUE), 'spikiness', 'shannon.entropy')
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
#' @importFrom stats dgeom dnbinom dpois reshape
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

	# Quality info
	if (is.null(model$qualityInfo$complexity)) { model$qualityInfo$complexity <- NA }
	if (is.null(model$qualityInfo$spikiness)) { model$qualityInfo$spikiness <- NA }
	if (is.null(model$qualityInfo$shannon.entropy)) { model$qualityInfo$shannon.entropy <- NA }
	if (is.null(model$qualityInfo$bhattacharyya)) { model$qualityInfo$bhattacharyya <- NA }
	quality.string <- paste0('reads = ',round(sum(model$bins$counts)/1e6,2),'M, complexity = ',round(model$qualityInfo$complexity),',  spikiness = ',round(model$qualityInfo$spikiness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2),',  bhattacharyya = ',round(model$qualityInfo$bhattacharyya,2), ', num.segments = ',length(model$segments))

	# Plot the histogram
	ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count") + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
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
			distributions[[length(distributions)+1]] <- weights[istate] * stats::dgeom(x, model$distributions[istate,'prob'])
		} else if (model$distributions[istate,'type']=='dnbinom') {
			# negative binomials
			distributions[[length(distributions)+1]] <- weights[istate] * stats::dnbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])
		} else if (model$distributions[istate,'type']=='dpois') {
			# poissons
			distributions[[length(distributions)+1]] <- weights[istate] * stats::dpois(x, model$distributions[istate,'lambda'])
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
	distributions <- stats::reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="state", times=c(c.state.labels,"total"))
	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions)
	} else {
		ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions[distributions$state==state,])
	}
	
	# Make legend and colors correct
	lmeans <- round(model$distributions[,'mu'], 2)
	lvars <- round(model$distributions[,'variance'], 2)
	lweights <- round(model$weights, 2)
	legend <- paste0(c.state.labels, ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
	legend <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=', round(var(counts),2), ', weight(data)=1'))
	ggplt <- ggplt + scale_color_manual(breaks=c(c.state.labels, 'total'), values=stateColors(c(c.state.labels,'total')), labels=legend)

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
#' @param plot.SCE Logical indicating whether SCE events should be plotted.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotKaryogram <- function(model, both.strands=FALSE, plot.SCE=FALSE, file=NULL) {

	if (class(model)=='GRanges') {
		binned.data <- model
		model <- list()
		model$ID <- ''
		model$bins <- binned.data
		model$qualityInfo <- list(shannon.entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity.preseqR'), bhattacharyya=NA)
		plot.karyogram(model, both.strands=both.strands, file=file)
	} else if (class(model)==class.univariate.hmm) {
		plot.karyogram(model, both.strands=both.strands, file=file)
	} else if (class(model)==class.bivariate.hmm) {
		plot.karyogram(model, both.strands=both.strands, plot.SCE=plot.SCE, file=file)
	}

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.karyogram <- function(model, both.strands=FALSE, plot.SCE=FALSE, file=NULL) {
	
	## Check user input
	if (is.null(model$sce) & plot.SCE) {
		warning("Cannot plot SCE coordinates. Please run 'getSCEcoordinates' first.")
		plot.SCE <- FALSE
	}

	## Convert to GRanges
	gr <- model$bins
	grl <- split(gr, seqnames(gr))

	## Get some variables
	fs.x <- 13
	num.chroms <- length(seqlevels(gr))
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

	# Quality info
	if (is.null(model$qualityInfo$complexity)) { model$qualityInfo$complexity <- NA }
	if (is.null(model$qualityInfo$spikiness)) { model$qualityInfo$spikiness <- NA }
	if (is.null(model$qualityInfo$shannon.entropy)) { model$qualityInfo$shannon.entropy <- NA }
	if (is.null(model$qualityInfo$bhattacharyya)) { model$qualityInfo$bhattacharyya <- NA }
	quality.string <- paste0('reads = ',round(sum(model$bins$counts)/1e6,2),'M, complexity = ',round(model$qualityInfo$complexity),',  spikiness = ',round(model$qualityInfo$spikiness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2),',  bhattacharyya = ',round(model$qualityInfo$bhattacharyya,2), ', num.segments = ',length(model$segments))

	## Get SCE coordinates
	if (plot.SCE) {
		scecoords <- model$sce
		# Set to midpoint
		start(scecoords) <- (start(scecoords)+end(scecoords))/2
		end(scecoords) <- start(scecoords)
	}

	## Theme for plotting chromosomes
	empty_theme <- theme(axis.line=element_blank(),
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		axis.title.x=element_text(size=fs.x),
		axis.title.y=element_blank(),
		legend.position="none",
		panel.background=element_blank(),
		panel.border=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		plot.background=element_blank())

	## Go through chromosomes and plot
	ggplts <- list()
	for (chrom in seqlevels(gr)) {
		i1 <- which(chrom==seqlevels(gr))

		# Plot the read counts
		dfplot <- as.data.frame(grl[[i1]])
		# Transform coordinates to match p-arm on top
		dfplot$start <- (-dfplot$start + seqlengths(gr)[i1])
		dfplot$end <- (-dfplot$end + seqlengths(gr)[i1])
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

		## Read counts
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
			ggplt <- ggplt + scale_color_manual(values=stateColors(levels(dfplot$state)), drop=FALSE)	# do not drop levels if not present
		} else {
			if (both.strands) {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts'), size=0.2, col=strandColors('+'))	# read count
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts'), size=0.2, col=strandColors('-'))	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')	# outliers
				ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')	# outliers
			} else {
				ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts'), size=0.2, col='gray20')	# read count
				ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts'), size=2, shape=21, col='gray20')	# outliers
			}
		}
		## Chromosome backbone
		if (both.strands) {
			ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=0, xmax=seqlengths(gr)[i1], col='white', fill='gray20')	# chromosome backbone as simple rectangle
		} else {
			ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, xmax=seqlengths(gr)[i1], col='white', fill='gray20')	# chromosome backbone as simple rectangle
		}
		if (plot.SCE) {
			dfsce <- as.data.frame(scecoords[seqnames(scecoords)==names(grl)[i1]])
			# Transform coordinates to match p-arm on top
			dfsce$start <- (-dfsce$start + seqlengths(gr)[i1])
			dfsce$end <- (-dfsce$end + seqlengths(gr)[i1])
			if (nrow(dfsce)>0) {
				ggplt <- ggplt + geom_segment(data=dfsce, aes(x=start, xend=start), y=-custom.xlim, yend=-0.5*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
			}
		}
		ggplt <- ggplt + empty_theme	# no axes whatsoever
		ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1]))	# chromosome names
		if (both.strands) {
			ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-custom.xlim,custom.xlim))	# set x- and y-limits
		} else {
			ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-0.6*custom.xlim,custom.xlim))	# set x- and y-limits
		}
		ggplts[[chrom]] <- ggplt
		
	}
	
	## Combine in one canvas
	fs.title <- 20
	nrows <- 2	# rows for plotting chromosomes
	nrows.text <- 2	# additional row for displaying ID and qualityInfo
	nrows.total <- nrows + nrows.text
	ncols <- ceiling(num.chroms/nrows)
	plotlist <- c(rep(list(NULL),ncols), ggplts, rep(list(NULL),ncols))
	cowplt <- plot_grid(plotlist=plotlist, nrow=nrows.total, rel_heights=c(2,21,21,2))
	cowplt <- cowplt + cowplot::draw_label(model$ID, x=0.5, y=0.99, vjust=1, hjust=0.5, size=fs.title)
	cowplt <- cowplt + cowplot::draw_label(quality.string, x=0.5, y=0.01, vjust=0, hjust=0.5, size=fs.x)
	if (!is.null(file)) {
		ggsave(file, cowplt, width=ncols*1.4, height=nrows*4.6)
	} else {
		return(cowplt)
	}

}


# =================================================================
# Plot a heatmap of chromosome state for multiple samples
# =================================================================
#' Plot aneuploidy state
#'
#' Plot a heatmap of aneuploidy state for multiple samples. Samples can be clustered and the output can be returned as data.frame.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmms}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity in their CNV-state.
#' @param as.data.frame If \code{TRUE}, instead of a plot, a data.frame with the aneuploidy state for each sample will be returned.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or a data.frame, depending on option \code{as.data.frame}.
#' @author Aaron Taudt
#' @importFrom stats aggregate dist hclust
#' @export
#'@examples
#'## Get results from a small-cell-lung-cancer
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'## Plot the ploidy state per chromosome
#'heatmapAneuploidies(files, cluster=FALSE)
#'## Return the ploidy state as data.frame
#'df <- heatmapAneuploidies(files, cluster=FALSE, as.data.frame=TRUE)
#'head(df)
#'
heatmapAneuploidies <- function(hmms, ylabels=NULL, cluster=TRUE, as.data.frame=FALSE) {

	## Check user input
	if (!is.null(ylabels)) {
		if (length(ylabels) != length(hmms)) {
			stop("length(ylabels) must equal length(hmms)")
		}
	}

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)
	levels.state <- unique(unlist(lapply(hmms, function(hmm) { levels(hmm$bins$state) })))
	
	## Assign new IDs
	if (!is.null(ylabels)) {
		for (i1 in 1:length(hmms)) {
			hmms[[i1]]$ID <- ylabels[i1]
		}
	}

	## Transform to GRanges in reduced representation
	grlred <- GRangesList()
	for (hmm in hmms) {
		if (!is.null(hmm$segments)) {
			grlred[[hmm$ID]] <- hmm$segments
		}
	}
	
	## Find the most frequent state (mfs) for each chromosome and sample
	ptm <- startTimedMessage("finding most frequent state for each sample and chromosome ...")
	grl.per.chrom <- lapply(grlred, function(x) { split(x, seqnames(x)) })
	mfs.samples <- list()
	for (i1 in 1:length(grlred)) {
		mfs.samples[[names(grlred)[i1]]] <- lapply(grl.per.chrom[[i1]], function(x) {
      if (length(x)>0) {
        tab <- stats::aggregate(width(x), by=list(state=x$state), FUN="sum")
        tab$state[which.max(tab$x)]
      } else {
        "0-somy"
      }
      })
		attr(mfs.samples[[names(grlred)[i1]]], "varname") <- 'chromosome'
	}
	attr(mfs.samples, "varname") <- 'sample'
	stopTimedMessage(ptm)

	## Transform to data.frame
	# Long format
	df <- reshape2::melt(mfs.samples, value.name='state')
	df$state <- factor(df$state, levels=levels.state)
	df$sample <- factor(df$sample, levels=unique(df$sample))
	df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
	# Wide format
	df.wide <- reshape2::dcast(df, sample ~ chromosome, value.var='state', factorsAsStrings=FALSE)
	# Correct strings to factors
	for (col in 2:ncol(df.wide)) {
		df.wide[,col] <- factor(df.wide[,col], levels=levels.state)
	}

	## Cluster the samples by chromosome state
	if (cluster) {
		# Cluster
		hc <- stats::hclust(stats::dist(data.matrix(df.wide[-1])))
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
			df.table[,i1] <- initializeStates(levels.state)$multiplicity[df.wide[,i1]]
		}
		return(df.table)
	} else {
		ggplt <- ggplot(df) + geom_tile(aes_string(x='chromosome', y='sample', fill='state'), col='black') + theme_bw() + scale_fill_manual(values=stateColors(levels(df$state)))
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
#' @param hmms A list of \code{\link{aneuHMM}} objects or files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmms}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param classes A character vector with the classification of the elements on the y-axis. The vector must have the same length as \code{hmms}. If specified the clustering algorithm will try to display similar categories together in the dendrogram.
#' @param classes.color A (named) vector with colors that are used to distinguish \code{classes}. Names must correspond to the unique elements in \code{classes}.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param plot.SCE Logical indicating whether SCE events should be plotted.
#' @param hotspots A \code{\link{GRanges}} object with coordinates of genomic hotspots (see \code{\link{hotspotter}}).
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @importFrom stats as.dendrogram
#' @export
#'@examples
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Plot a clustered heatmap
#'classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#'labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#'heatmapGenomewide(c(lung.files, liver.files), ylabels=labels, classes=classes,
#'                  classes.color=c('blue','red'))
#'
heatmapGenomewide <- function(hmms, ylabels=NULL, classes=NULL, classes.color=NULL, file=NULL, cluster=TRUE, plot.SCE=TRUE, hotspots=NULL) {

	## Check user input
	if (!is.null(ylabels)) {
		if (length(ylabels) != length(hmms)) {
			stop("length(ylabels) must equal length(hmms)")
		}
	}
	if (!is.null(classes)) {
		if (length(classes) != length(hmms)) {
			stop("length(classes) must equal length(hmms)")
		}
	}
	if (length(classes.color)!=length(unique(classes))) {
		stop("'classes.color' must have the same length as unique(classes)")
	}
	if (is.null(names(classes.color))) {
		names(classes.color) <- unique(classes)
	}
	if (!setequal(names(classes.color), unique(classes))) {
		stop("The names of 'classes.color' must be equal to the unique elements in 'classes'")
	}

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)

	## Dataframe with IDs, ylabels and classes
	data <- data.frame(ID=sapply(hmms,'[[','ID'))
	data$ylabel <- ylabels
	data$class <- classes

	## Assign new IDs
	if (!is.null(ylabels)) {
		for (i1 in 1:length(hmms)) {
			hmms[[i1]]$ID <- as.character(ylabels[i1])
		}
	}

	## Get segments and SCE coordinates
	temp <- getSegments(hmms, cluster=cluster, classes=classes)
	grlred <- temp$segments
	hc <- temp$clustering
	if (cluster) {
		hmms <- hmms[hc$order]
		data <- data[hc$order,]
	}
	if (plot.SCE) {
		sce <- lapply(hmms,'[[','sce')
		sce <- sce[!unlist(lapply(sce, is.null))]
		sce <- sce[lapply(sce, length)!=0]
		if (length(sce)==0) {
			plot.SCE <- FALSE
		}
	}

	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
	ptm <- startTimedMessage("Transforming coordinates ...")
	grlred <- endoapply(grlred, transCoord)
	if (plot.SCE) {
		sce <- endoapply(sce, transCoord)
	}
	stopTimedMessage(ptm)

	## Data.frame for plotting
	ptm <- startTimedMessage("Making the plot ...")
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
	pltlist <- list()
	widths <- vector()

	## Prepare the plot
	ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=stateColors(levels(df$state))) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20), axis.line=element_blank())
	ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
	if (plot.SCE) {
		ggplt <- ggplt + geom_linerange(data=df.sce, mapping=aes_string(x='sample', ymin='start', ymax='end'), size=2) + ylab('')
	}
	if (!is.null(hotspots)) {
		df.hot <- as.data.frame(transCoord(hotspots))
		df.hot$xmin <- 0
		df.hot$xmax <- length(unique(df$sample))+1
		ggplt <- ggplt + geom_rect(data=df.hot, mapping=aes_string(xmin='xmin', xmax='xmax', ymin='start.genome', ymax='end.genome', alpha='num.events'), fill='hotpink4') + scale_alpha_continuous(name='SCE events', range=c(0.4,0.8))
	}
	width.heatmap <- sum(as.numeric(seqlengths(hmms[[1]]$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
	height <- length(hmms) * 0.5
	pltlist[['heatmap']] <- ggplt
	widths['heatmap'] <- width.heatmap
	## Make the classification bar
	if (!is.null(classes)) {
		width.classes <- 5
		data$y <- 1:nrow(data)
		ggclass <- ggplot(data) + geom_tile(aes_string(x=1, y='y', fill='class')) + guides(fill=FALSE) + theme(axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) + coord_cartesian(ylim=c(1.5,nrow(data)-0.5))
		if (!is.null(classes.color)) {
			ggclass <- ggclass + scale_fill_manual(breaks=names(classes.color), values=classes.color)
		}
		pltlist[['classbar']] <- ggclass
		widths['classbar'] <- width.classes
	}
	## Prepare the dendrogram
	if (!is.null(hc)) {
		dhc <- stats::as.dendrogram(hc)
		ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
		ggdndr <- ggplot(ggdendro::segment(ddata)) + geom_segment(aes_string(x='x', y='y', xend='xend', yend='yend')) + coord_flip(xlim=c(1.5,nrow(ddata$labels)-0.5)) + scale_y_reverse(expand=c(0,0)) + ggdendro::theme_dendro()
		width.dendro <- 20
		pltlist[['dendro']] <- ggdndr
		widths['dendro'] <- width.dendro
	}
	cowplt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths=rev(widths))
	stopTimedMessage(ptm)

	## Plot to file
	if (!is.null(file)) {
		ptm <- startTimedMessage("plotting to file ",file," ...")
		ggsave(file, cowplt, width=sum(widths)/2.54, height=height/2.54, limitsize=FALSE)
		stopTimedMessage(ptm)
	} else {
		return(cowplt)
	}

}


# =================================================================
# Plot profile with read counts
# =================================================================
#' Read count and CNV profile
#'
#' Plot a profile with read counts and CNV-state from a \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#'
#' @param model A \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#' @param file A PDF file where the plot will be saved.
#' @param plot.SCE Logical indicating whether SCE events should be plotted.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotProfile <- function(model, both.strands=FALSE, plot.SCE=TRUE, file=NULL) {

	if (class(model)=='GRanges') {
		binned.data <- model
		model <- list()
		model$ID <- ''
		model$bins <- binned.data
		model$qualityInfo <- list(shannon.entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity.preseqR'))
		plot.profile(model, both.strands=both.strands, plot.SCE=FALSE, file=file)
	} else if (class(model)==class.univariate.hmm) {
		plot.profile(model, both.strands=FALSE, plot.SCE=FALSE, file=file)
	} else if (class(model)==class.bivariate.hmm) {
		plot.profile(model, both.strands=both.strands, plot.SCE=plot.SCE, file=file)
	}

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.profile <- function(model, both.strands=FALSE, plot.SCE=TRUE, file=NULL) {
	
	## Convert to GRanges
	bins <- model$bins
	## Get SCE coordinates
	if (is.null(model$sce) & plot.SCE) {
		warning("Cannot plot SCE coordinates. Please run 'getSCEcoordinates' first.")
		plot.SCE <- FALSE
	}
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
# 	custom.xlim1 <- as.numeric(names(tab)[which.max(tab)]) # maximum value of read distribution
# 	custom.xlim2 <- as.integer(mean(bins$counts, trim=0.05)) # mean number of counts
# 	custom.xlim <- max(custom.xlim1, custom.xlim2, na.rm=TRUE) * 2.7
	custom.xlim <- get_rightxlim(bins$counts)
	if (both.strands) {
		custom.xlim <- custom.xlim / 1
	}
	if (length(custom.xlim)==0 | is.na(custom.xlim)) {
		custom.xlim <- 1
	}

	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
	ptm <- startTimedMessage("transforming coordinates ...")
	bins <- transCoord(bins)
	if (plot.SCE) {
		scecoords <- transCoord(scecoords)
	}
	stopTimedMessage(ptm)

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
		axis.ticks=element_blank(),
		axis.title.x=element_blank(),
		panel.background=element_blank(),
		panel.border=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		plot.background=element_blank())
	ggplt <- ggplot(dfplot, aes_string(x='start.genome', y='counts'))	# data
	if (both.strands) {
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='pcounts'), position=position_jitter(width=0, height=0))	# read count
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='mcounts'), position=position_jitter(width=0, height=0))	# read count
# 		ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='pcounts'))	# read count
# 		ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='mcounts'))	# read count
	} else {
		ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='counts'), position=position_jitter(width=0, height=0))	# read count
# 		ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='counts'))	# read count
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
	
	ggplt <- ggplt + scale_color_manual(values=stateColors(levels(dfplot$state)), drop=FALSE)	# do not drop levels if not present
	if (plot.SCE) {
		dfsce <- as.data.frame(scecoords)
		if (nrow(dfsce)>0) {
			ggplt <- ggplt + geom_segment(data=dfsce, aes_string(x='start.genome', xend='start.genome'), y=-1.5*custom.xlim, yend=-1.3*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
		}
	}
	ggplt <- ggplt + empty_theme	# no axes whatsoever
	if (both.strands) {
		ggplt <- ggplt + coord_cartesian(ylim=c(-1.5*custom.xlim,custom.xlim))	# set x- and y-limits
	} else {
		ggplt <- ggplt + coord_cartesian(ylim=c(0,custom.xlim))	# set x- and y-limits
	}
	# Get midpoints of each chromosome for xticks
	ggplt <- ggplt + scale_x_continuous(breaks=seqlengths(model$bins)/2+cum.seqlengths.0[as.character(seqlevels(model$bins))], labels=seqlevels(model$bins))
	# Quality info
	if (is.null(model$qualityInfo$complexity)) { model$qualityInfo$complexity <- NA }
	if (is.null(model$qualityInfo$spikiness)) { model$qualityInfo$spikiness <- NA }
	if (is.null(model$qualityInfo$shannon.entropy)) { model$qualityInfo$shannon.entropy <- NA }
	if (is.null(model$qualityInfo$bhattacharyya)) { model$qualityInfo$bhattacharyya <- NA }
	quality.string <- paste0('reads = ',round(sum(model$bins$counts)/1e6,2),'M, complexity = ',round(model$qualityInfo$complexity),',  spikiness = ',round(model$qualityInfo$spikiness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2),',  bhattacharyya = ',round(model$qualityInfo$bhattacharyya,2), ', num.segments = ',length(model$segments))
	ggplt <- ggplt + ylab('read count') + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
		
	if (!is.null(file)) {
		ggsave(file, ggplt, width=20, height=5)
	} else {
		return(ggplt)
	}
}

