

#' @import ggplot2
#' @importFrom cowplot plot_grid draw_label
#' @import reshape2
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
    state.colors <- c("zero-inflation"="gray90", "0-somy"="gray90","1-somy"="darkorchid3","2-somy"="springgreen2","3-somy"="red3","4-somy"="gold2","5-somy"="navy","6-somy"="lemonchiffon","7-somy"="dodgerblue","8-somy"="chartreuse4","9-somy"="lightcoral","10-somy"="aquamarine2","total"="black")
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

#' @describeIn colors Colors that are used for breakpoint types.
#' @param breaktypes A character vector with breakpoint types whose color should be returned. Any combination of \code{c('CNB','SCE','CNB+SCE','other')}.
#' @export
#'@examples
#'## Make a nice pie chart with the AneuFinder breakpoint-type color scheme
#'breakpointcolors <- breakpointColors()
#'pie(rep(1,length(breakpointcolors)), labels=names(breakpointcolors), col=breakpointcolors)
#'
breakpointColors <- function(breaktypes=c('CNB','SCE','CNB+SCE','other')) {
    break.colors <- c('CNB'="dodgerblue4", 'SCE'="tomato3", 'CNB+SCE'="orchid4", 'other'='gray30')
    breaks.with.color <- intersect(breaktypes, names(break.colors))
    cols <- rep('gray30', length(breaktypes))
    names(cols) <- breaktypes
    cols[breaks.with.color] <- break.colors[breaks.with.color]
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
#' @param x A \code{\link{GRanges-class}} object with binned read counts.
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
        plotHistogram(x, ...)
    } else if (type == 'profile' | type==1) {
        plotProfile(x, ...)
    }
}

#' Plotting function for binned read counts (list)
#'
#' Make plots for binned read counts (list) from \code{\link{binned.data}}.
#'
#' @param x A \code{\link{GRangesList}} object with binned read counts.
#' @param type Type of the plot, one of \code{c('profile', 'histogram', 'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with read counts.}
#'   \item{\code{histogram}}{A histogram of read counts.}
#'   \item{\code{profile}}{An profile with read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRangesList
#' @export
plot.GRangesList <- function(x, type='profile', ...) {
    if (type == 'karyogram' | type==3) {
        plotKaryogram(x, ...)
    } else if (type == 'histogram' | type==2) {
        plotHistogram(x, ...)
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
        plotHistogram(x, ...)
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
#     rightxlim1 <- median(counts[counts>0])*7
#     tab <- table(counts)
#     tab <- tab[names(tab)!='0']
#     breaks <- as.numeric(names(tab))
#     rightxlim2 <- breaks[tab<=5 & breaks>median(counts)*2][1]
#     rightxlim <- min(rightxlim1,rightxlim2, na.rm=TRUE)
    rightxlim <- stats::quantile(counts, 0.999)
    if (length(rightxlim)==0 | is.na(rightxlim) | is.infinite(rightxlim)) {
        rightxlim <- 1
    }
    return(rightxlim)
}

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} with two additional metadata columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
    cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(gr)
    gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    return(gr)
}

# =================================================================
# Plot a read histogram with univariate fits for a bivariate HMM
# =================================================================
plotBivariateHistograms <- function(bihmm) {

    ## Stack the two strands (bins)
    binned.data <- bihmm$bins
    
    binned.data.minus <- binned.data
    strand(binned.data.minus) <- '-'
    mcols(binned.data.minus) <- NULL
    binned.data.minus$state <- binned.data$mstate
    binned.data.minus$copy.number <- binned.data$mcopy.number
    
    binned.data.plus <- binned.data
    strand(binned.data.plus) <- '+'
    mcols(binned.data.plus) <- NULL
    binned.data.plus$state <- binned.data$pstate
    binned.data.plus$copy.number <- binned.data$pcopy.number
    
    binned.data.stacked <- c(binned.data.minus, binned.data.plus)
    
    ## Attributes
    mask.attributes <- c('complexity', 'spikiness', 'entropy')
    attributes(binned.data.stacked)[mask.attributes] <- attributes(binned.data)[mask.attributes]
    
    ## Stack the two strands (bincounts)
    bincounts <- bihmm$bincounts[[1]]
    
    bincounts.minus <- bincounts
    strand(bincounts.minus) <- '-'
    mcols(bincounts.minus) <- NULL
    bincounts.minus$counts <- bincounts$mcounts
    
    bincounts.plus <- bincounts
    strand(bincounts.plus) <- '+'
    mcols(bincounts.plus) <- NULL
    bincounts.plus$counts <- bincounts$pcounts
    
    bincounts.stacked <- c(bincounts.minus, bincounts.plus)

    ## Make fake uni.hmm and plot
    strand <- 'minus'
    uni.hmm <- list()
    uni.hmm$ID <- bihmm$ID
    uni.hmm$bins <- binned.data.stacked
    uni.hmm$bincounts <- GRangesList(bincounts.stacked)
    uni.hmm$segments <- bihmm$segments
    uni.hmm$weights <- bihmm$univariateParams$weights
    uni.hmm$distributions <- bihmm$distributions[[strand]]
    uni.hmm$qualityInfo <- bihmm$qualityInfo
    class(uni.hmm) <- "aneuHMM"
    ggplts <- plotHistogram(uni.hmm)

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
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom stats dgeom dnbinom dbinom dpois reshape
plotHistogram <- function(model) {

    model <- suppressMessages( loadFromFiles(model, check.class=c('GRanges', 'GRangesList', "aneuHMM"))[[1]] )
    if (class(model) == 'GRanges') {
        bins <- model
        bincounts <- model
    } else if (is(model, "GRangesList")) {
        bins <- model[[1]]
        bincounts <- model[[1]]
    } else if (class(model) == "aneuHMM") {
        bins <- model$bins
        if (!is.null(model$bins$counts)) {
            bincounts <- model$bins
        } else if (!is.null(model$bincounts[[1]]$counts)) {
            bincounts <- model$bincounts[[1]]
        }
    }
    counts <- bincounts$counts

      states <- bins$state
        if (!is.null(model$weights)) {
              weights <- model$weights
        }

      # Find the x limits
      breaks <- max(counts)
      if (max(counts)==0) { breaks <- 1 }
      rightxlim <- get_rightxlim(counts)

      # Quality info
      qualityInfo <- getQC(model)
      quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))

      # Plot the histogram
      ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count") + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
      if (is.null(model$weights)) {
            return(ggplt)
      }

      if (!is.null(model$bins$state)) {
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
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dpois(x, model$distributions[istate,'prob'])
                } else if (model$distributions[istate,'type']=='dbinom') {
                      # binomials
                      s <- model$distributions[istate,'size']
                      p <- model$distributions[istate,'prob']
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dbinom(x, round(s), p)
                }
          }
          distributions <- as.data.frame(distributions)
          names(distributions) <- c("x",c.state.labels)
          # Total
          distributions$total <- apply(distributions[-1], 1, sum)
      
          # Reshape the data.frame for plotting with ggplot
          distributions <- stats::reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="state", times=c(c.state.labels,"total"))
          ### Plot the distributions
            ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions)
          
          # Make legend and colors correct
          lmeans <- round(model$distributions[,'mu'], 2)
          lvars <- round(model$distributions[,'variance'], 2)
          lweights <- round(model$weights, 2)
          legend <- paste0(c.state.labels, ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
          legend <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=', round(var(counts),2), ', weight(data)=1'))
          ggplt <- ggplt + scale_color_manual(breaks=c(c.state.labels, 'total'), values=stateColors(c(c.state.labels,'total')), labels=legend)
      }
  
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
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotKaryogram <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL) {

    if (class(model)=='GRanges') {
        binned.data <- model
        model <- list()
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'), bhattacharyya=NA)
        plot.karyogram(model, both.strands=both.strands, file=file)
    } else if (class(model)=="aneuHMM") {
        plot.karyogram(model, both.strands=both.strands, file=file)
    } else if (class(model)=="aneuBiHMM") {
        plot.karyogram(model, both.strands=both.strands, plot.breakpoints=plot.breakpoints, file=file)
    }

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.karyogram <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL) {
    
    ## Check user input
    if (is.null(model$breakpoints) & plot.breakpoints) {
        warning("Cannot find breakpoint coordinates. Please run 'getBreakpoints' first.")
        plot.breakpoints <- FALSE
    }

    ## Convert to GRanges
  if (!is.null(model$bins$counts)) {
    bins <- model$bins
  } else if (!is.null(model$bincounts[[1]]$counts)) {
    bins <- model$bincounts[[1]]
    ind <- findOverlaps(bins, model$bins, select='first')
    bins$state <- model$bins$state[ind]
    bins$mstate <- model$bins$mstate[ind]
    bins$pstate <- model$bins$pstate[ind]
  }
    bins.split <- split(bins, seqnames(bins))
    bins.split <- bins.split[lengths(bins.split) > 0]

    ## Get some variables
    fs.x <- 13
    maxseqlength <- max(seqlengths(bins))
    tab <- table(bins$counts)
    tab <- tab[names(tab)!='0']
    if (both.strands) {
      custom.xlim <- get_rightxlim(c(bins$mcounts, bins$pcounts))
    } else {
      custom.xlim <- get_rightxlim(bins$counts)
    }

    # Quality info
    qualityInfo <- getQC(model)
    quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))

    ## Get breakpoint coordinates
    if (plot.breakpoints) {
        bp.coords <- model$breakpoints
        # Set to midpoint
        start(bp.coords) <- (start(bp.coords)+end(bp.coords))/2
        end(bp.coords) <- start(bp.coords)
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
    for (i1 in 1:length(bins.split)) {
      chrom <- names(bins.split)[i1]

        # Plot the read counts
        dfplot <- as.data.frame(bins.split[[i1]])
        # Transform coordinates to match p-arm on top
        dfplot$start <- (-dfplot$start + seqlengths(bins)[chrom])
        dfplot$end <- (-dfplot$end + seqlengths(bins)[chrom])
        # Set values too big for plotting to limit
        dfplot$counts[dfplot$counts>=custom.xlim] <- custom.xlim
        dfplot.points <- dfplot[dfplot$counts>=custom.xlim,]
        dfplot.points$counts <- rep(custom.xlim, nrow(dfplot.points))

        if (both.strands) {
            dfplot$mcounts <- - dfplot$mcounts    # negative minus counts
            dfplot$pcounts[dfplot$pcounts>=custom.xlim] <- custom.xlim
            dfplot$mcounts[dfplot$mcounts<=-custom.xlim] <- -custom.xlim
            dfplot.points.plus <- dfplot[dfplot$pcounts>=custom.xlim,]
            dfplot.points.plus$counts <- rep(custom.xlim, nrow(dfplot.points.plus))
            dfplot.points.minus <- dfplot[dfplot$mcounts<=-custom.xlim,]
            dfplot.points.minus$counts <- rep(-custom.xlim, nrow(dfplot.points.minus))
        }

        ## Read counts
        ggplt <- ggplot(dfplot, aes_string(x='start', y='counts'))    # data
        if (!is.null(bins.split[[i1]]$state)) {
            if (both.strands) {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts', col='pstate'), size=0.2)    # read count
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts', col='mstate'), size=0.2)    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts', col='pstate'), size=5, shape=21)    # outliers
                ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts', col='mstate'), size=5, shape=21)    # outliers
            } else {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts', col='state'), size=0.2)    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts', col='state'), size=2, shape=21)    # outliers
            }
          statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
            ggplt <- ggplt + scale_color_manual(values=stateColors(statelevels), drop=FALSE)    # do not drop levels if not present
        } else {
            if (both.strands) {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts'), size=0.2, col=strandColors('+'))    # read count
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts'), size=0.2, col=strandColors('-'))    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')    # outliers
                ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')    # outliers
            } else {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts'), size=0.2, col='gray20')    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts'), size=2, shape=21, col='gray20')    # outliers
            }
        }
        ## Chromosome backbone
        if (both.strands) {
            ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=0, xmax=seqlengths(bins)[chrom], col='white', fill='gray20')    # chromosome backbone as simple rectangle
        } else {
            ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, xmax=seqlengths(bins)[chrom], col='white', fill='gray20')    # chromosome backbone as simple rectangle
        }
        if (plot.breakpoints) {
            df.bp <- as.data.frame(bp.coords[seqnames(bp.coords)==names(bins.split)[i1]])
            # Transform coordinates to match p-arm on top
            df.bp$start <- (-df.bp$start + seqlengths(bins)[chrom])
            df.bp$end <- (-df.bp$end + seqlengths(bins)[chrom])
            if (nrow(df.bp)>0) {
            statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
              suppressMessages( ggplt <- ggplt + scale_color_manual(values=c(breakpointColors(), stateColors(statelevels)), drop=FALSE) )
                ggplt <- ggplt + geom_segment(data=df.bp, aes_string(x='start', xend='start', color='type'), y=-custom.xlim, yend=-0.5*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'), alpha=0.5)
            }
        }
        ggplt <- ggplt + empty_theme    # no axes whatsoever
        ggplt <- ggplt + ylab(names(bins.split)[i1])    # chromosome names
        if (both.strands) {
            ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-custom.xlim,custom.xlim))    # set x- and y-limits
        } else {
            ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-0.6*custom.xlim,custom.xlim))    # set x- and y-limits
        }
        ggplts[[i1]] <- ggplt
        
    }
    names(ggplts) <- names(bins.split)
    
    ## Combine in one canvas
    fs.title <- 20
    nrows <- 2    # rows for plotting chromosomes
    nrows.text <- 2    # additional row for displaying ID and qualityInfo
    nrows.total <- nrows + nrows.text
    ncols <- ceiling(length(bins.split)/nrows)
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
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
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
  if (length(hmms) == 1 & cluster==TRUE) {
    cluster <- FALSE
    warning("Cannot do clustering because only one object was given.")
  }

    ## Load the files
    hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
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
        mfs.samples[[names(grlred)[i1]]] <- list()
        for (i2 in 1:length(grl.per.chrom[[i1]])) {
          x <- grl.per.chrom[[i1]][[i2]]
      if (length(x)>0) {
        tab <- stats::aggregate(width(x), by=list(state=x$state), FUN="sum")
            mfs.samples[[names(grlred)[i1]]][[i2]] <- tab$state[which.max(tab$x)]
      } else {
            mfs.samples[[names(grlred)[i1]]][[i2]] <- "0-somy"
      }
        }
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
        ## Reorder state levels for the legend
        df$state <- factor(df$state, levels=names(sort(initializeStates(levels(df$state))$multiplicity)))
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
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the same length as \code{hmms}. If \code{NULL} the IDs from the \code{\link{aneuHMM}} objects will be used.
#' @param classes A character vector with the classification of the elements on the y-axis. The vector must have the same length as \code{hmms}.
#' @param reorder.by.class If \code{TRUE}, the dendrogram will be reordered to display similar classes next to each other.
#' @param classes.color A (named) vector with colors that are used to distinguish \code{classes}. Names must correspond to the unique elements in \code{classes}.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the samples should be clustered by similarity in their CNV-state.
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @param hotspots A \code{\link{GRanges-class}} object with coordinates of genomic hotspots (see \code{\link{hotspotter}}).
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will be excluded from the computation of the clustering. This can be useful to exclude regions with artifacts.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @importFrom stats as.dendrogram
#' @importFrom ggdendro dendro_data theme_dendro
#' @importFrom S4Vectors endoapply
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
heatmapGenomewide <- function(hmms, ylabels=NULL, classes=NULL, reorder.by.class=TRUE, classes.color=NULL, file=NULL, cluster=TRUE, plot.breakpoints=FALSE, hotspots=NULL, exclude.regions=NULL) {

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
  if (length(hmms) == 1 & cluster==TRUE) {
    cluster <- FALSE
    warning("Cannot do clustering because only one object was given.")
  }

    ## Load the files
    hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))

    ## Dataframe with IDs, ylabels and classes
    class.data <- data.frame(ID=sapply(hmms,'[[','ID'))
    class.data$ID <- factor(class.data$ID, levels=class.data$ID)
    if (is.null(ylabels)) {
      class.data$ylabel <- as.character(class.data$ID)
    } else {
      class.data$ylabel <- as.character(ylabels)
    }
    class.data$class <- classes
    
    ## Mapping to match ID with ylabel
    mapping <- class.data$ylabel
    names(mapping) <- class.data$ID

    ## Cluster
    if (reorder.by.class) {
      cl <- clusterHMMs(hmms, cluster=cluster, classes=classes, exclude.regions = exclude.regions)
    } else {
      cl <- clusterHMMs(hmms, cluster=cluster, exclude.regions = exclude.regions)
    }
    hmms <- hmms[cl$IDorder]
    class.data <- class.data[cl$IDorder,]
    class.data$ID <- factor(class.data$ID, levels=class.data$ID)
    ## Extract segements
  segments.list <- GRangesList()
  for (i1 in 1:length(hmms)) {
      hmm <- hmms[[i1]]
      if (is.null(hmm$segments)) {
          segments.list[[hmm$ID]] <- GRanges()
      } else {
            segments.list[[hmm$ID]] <- hmm$segments
      }
  }
  ## Extract breakpoints    
    if (plot.breakpoints) {
      breakpoints <- GRangesList()
      for (i1 in 1:length(hmms)) {
        hmm <- hmms[[i1]]
          if (is.null(hmm$breakpoints)) {
              breakpoints[[hmm$ID]] <- GRanges()
          } else {
              breakpoints[[hmm$ID]] <- hmm$breakpoints
          }
      }
        if (length(breakpoints)==0) {
            plot.breakpoints <- FALSE
        }
    }

    ## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
    ptm <- startTimedMessage("Transforming coordinates ...")
    segments.list <- endoapply(segments.list, transCoord)
    if (plot.breakpoints) {
        breakpoints <- endoapply(breakpoints, transCoord)
    }
    stopTimedMessage(ptm)

    ## Data.frame for plotting
    ptm <- startTimedMessage("Making the plot ...")
    df <- list()
    for (i1 in 1:length(segments.list)) {
        df[[length(df)+1]] <- data.frame(start=segments.list[[i1]]$start.genome, end=segments.list[[i1]]$end.genome, seqnames=seqnames(segments.list[[i1]]), ID=names(segments.list)[i1], state=segments.list[[i1]]$state)
    }
    df <- do.call(rbind, df)
    df$ID <- factor(df$ID, levels=levels(class.data$ID))
    df$ylabel <- mapping[as.character(df$ID)]
    
    if (plot.breakpoints) {
        df.breakpoints <- list()
        for (i1 in 1:length(breakpoints)) {
          if (length(breakpoints[[i1]]) > 0) {
                df.breakpoints[[length(df.breakpoints)+1]] <- data.frame(start=breakpoints[[i1]]$start.genome, end=breakpoints[[i1]]$end.genome, seqnames=seqnames(breakpoints[[i1]]), ID=names(segments.list)[i1], mid=(breakpoints[[i1]]$start.genome + breakpoints[[i1]]$end.genome)/2)
          } else {
                df.breakpoints[[length(df.breakpoints)+1]] <- data.frame(start=numeric(), end=numeric(), seqnames=character(), ID=character(), mid=numeric())
          }
        }
        df.breakpoints <- do.call(rbind, df.breakpoints)
      df.breakpoints$ID <- factor(df.breakpoints$ID, levels=levels(class.data$ID))
      df.breakpoints$ylabel <- mapping[as.character(df.breakpoints$ID)]
    }
    
    # Chromosome lines
    cum.seqlengths <- cumsum(as.numeric(seqlengths(segments.list[[1]])))
    names(cum.seqlengths) <- seqlevels(segments.list[[1]])
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(segments.list[[1]])
    label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(segments.list[[1]]) )
    df.chroms <- data.frame(y=c(0,cum.seqlengths), x=1, xend=length(segments.list))

    ### Plot ###
    pltlist <- list()
    widths <- vector()

    ## Prepare the plot
    df$state <- factor(df$state, levels=names(sort(initializeStates(levels(df$state))$multiplicity)))
    df$x <- as.numeric(df$ID) # transform all x-coordiantes to numeric because factors and numerics get selected different margins
    ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='x', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + scale_x_continuous(name="sample", breaks=1:length(unique(df$ylabel)), labels=unique(df$ylabel))
    ggplt <- ggplt + scale_color_manual(values=stateColors(levels(df$state)))
    ggplt <- ggplt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20), axis.line=element_blank(), axis.title.x=element_blank())
    # ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
    ggplt <- ggplt + geom_segment(aes_string(x='x', xend='xend', y='y', yend='y'), data=df.chroms, col='black')
    ggplt <- ggplt + coord_flip()
    if (plot.breakpoints) {
        df.breakpoints$x <- as.numeric(df.breakpoints$ID)
        ggplt <- ggplt + geom_linerange(data=df.breakpoints, mapping=aes_string(x='x', ymin='start', ymax='end'), size=2) + ylab('') + geom_point(data=df.breakpoints, mapping=aes_string(x='x', y='mid'))
    }
    if (!is.null(hotspots)) {
      if (length(hotspots) > 0) {
          df.hot <- as.data.frame(transCoord(hotspots))
          df.hot$xmin <- 0
          df.hot$xmax <- length(class.data$ID)+1
          ggplt <- ggplt + geom_rect(data=df.hot, mapping=aes_string(xmin='xmin', xmax='xmax', ymin='start.genome', ymax='end.genome', alpha='num.events'), fill='hotpink4') + scale_alpha_continuous(name='breakpoints', range=c(0.4,0.8))
      }
    }
    width.heatmap <- sum(as.numeric(seqlengths(hmms[[1]]$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
    height <- max(length(hmms) * 0.5, 2)
    pltlist[['heatmap']] <- ggplt
    widths['heatmap'] <- width.heatmap
    ## Make the classification bar
    if (!is.null(classes)) {
        width.classes <- 5
      class.data$x <- as.numeric(class.data$ID)  # transform all x-coordiantes to numeric because factors and numerics get selected different margins
        ggclass <- ggplot(class.data) + geom_linerange(aes_string(ymin=0, ymax=1, x='x', col='class'), size=5) + guides(col=FALSE) + xlab("class")
      ggclass <- ggclass + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title.x=element_blank())
        ggclass <- ggclass + coord_flip()
        if (!is.null(classes.color)) {
            ggclass <- ggclass + scale_color_manual(breaks=names(classes.color), values=classes.color)
        }
        pltlist[['classbar']] <- ggclass
        widths['classbar'] <- width.classes
    }
    ## Prepare the dendrogram
    if (!is.null(cl$hclust)) {
        dhc <- stats::as.dendrogram(cl$hclust)
        ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
        ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x='x', xend='xend', y='y', yend='yend')) + scale_y_reverse()
        ggdndr <- ggdndr + coord_flip()
        ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
        width.dendro <- 20
        pltlist[['dendro']] <- ggdndr
        widths['dendro'] <- width.dendro
    }
    cowplt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths=rev(widths))
    stopTimedMessage(ptm)

    ## Plot to file
    if (!is.null(file)) {
        ptm <- startTimedMessage("Plotting to file ",file," ...")
        ggsave(file, cowplt, width=sum(widths), height=height, units='cm', limitsize=FALSE)
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
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @param normalize.counts An character giving the copy number state to which to normalize the counts, e.g. '1-somy', '2-somy' etc.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotProfile <- function(model, both.strands=FALSE, plot.breakpoints=FALSE, file=NULL, normalize.counts=NULL) {

    if (class(model)=='GRanges') {
        binned.data <- model
        model <- list()
        class(model) <- "aneuHMM"
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'))
        plot.profile(model, both.strands=both.strands, plot.breakpoints=FALSE, file=file)
    } else if (is(model, "GRangesList")) {
        binned.data <- model[[1]]
        model <- list()
        class(model) <- "aneuHMM"
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'))
        plot.profile(model, both.strands=both.strands, plot.breakpoints=FALSE, file=file)
    } else if (class(model)=="aneuHMM") {
        plot.profile(model, both.strands=FALSE, plot.breakpoints=FALSE, file=file, normalize.counts = normalize.counts)
    } else if (class(model)=="aneuBiHMM") {
        plot.profile(model, both.strands=both.strands, plot.breakpoints=plot.breakpoints, file=file, normalize.counts = normalize.counts)
    }

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.profile <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL, normalize.counts=NULL) {
    
    ## Convert to GRanges
    if (!is.null(model$bins$counts)) {
        bins <- model$bins
    } else if (!is.null(model$bincounts[[1]]$counts)) {
        bins <- model$bincounts[[1]]
    }
    ## Get breakpoint coordinates
    if (is.null(model$breakpoints) & plot.breakpoints) {
        warning("Cannot breakpoints coordinates. Please run 'getBreakpoints' first.")
        plot.breakpoints <- FALSE
    }
    if (plot.breakpoints) {
        bp.coords <- model$breakpoints
        # Set to midpoint
        start(bp.coords) <- (start(bp.coords)+end(bp.coords))/2
        end(bp.coords) <- start(bp.coords)
    }

    ## Get some variables
    num.chroms <- length(levels(seqnames(bins)))
    maxseqlength <- max(seqlengths(bins))
    tab <- table(bins$counts)
    tab <- tab[names(tab)!='0']
    if (both.strands) {
      custom.xlim <- get_rightxlim(c(bins$mcounts, bins$pcounts))
    } else {
      custom.xlim <- get_rightxlim(bins$counts)
    }

    ## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
    bins <- transCoord(bins)
    if (plot.breakpoints) {
        bp.coords <- transCoord(bp.coords)
    }

    # Plot the read counts
    dfplot <- as.data.frame(bins)
    # Set values too big for plotting to limit
    if (both.strands) {
        dfplot$mcounts <- - dfplot$mcounts    # negative minus counts
    }
    # Mean counts for CNV-state
    if (!is.null(model$segments$state)) {
        dfplot.seg <- as.data.frame(transCoord(model$segments))
        if (class(model)=="aneuHMM") {
            dfplot.seg$counts.CNV <- model$distributions[as.character(dfplot.seg$state),'mu']
        } else if (class(model)=="aneuBiHMM") {
            if (!is.null(model$distributions$both)) {
                dfplot.seg$counts.CNV <- model$distributions$both[as.character(dfplot.seg$state),'mu']
            } else {
                dfplot.seg$counts.CNV <- model$distributions$plus[as.character(dfplot.seg$state),'mu']
            }
            dfplot.seg$pcounts.CNV <- model$distributions$plus[as.character(dfplot.seg$pstate),'mu']
            dfplot.seg$mcounts.CNV <- -model$distributions$minus[as.character(dfplot.seg$mstate),'mu']
        }
    }
    # Normalize counts
    ylabstring <- 'read count'
    if (!is.null(normalize.counts)) {
        colmask <- grepl('counts', names(dfplot))
        if (class(model)=="aneuHMM") {
            nfactor <- model$distributions[normalize.counts, 'mu']
        } else if (class(model)=="aneuBiHMM") {
            nfactor <- model$distributions$plus[normalize.counts, 'mu']
        }
        dfplot[,colmask] <- dfplot[,colmask] / nfactor
        colmask <- grepl('counts', names(dfplot.seg))
        dfplot.seg[,colmask] <- dfplot.seg[,colmask] / nfactor
        custom.xlim <- custom.xlim / nfactor
        ylabstring <- 'normalized read count'
    }

    empty_theme <- theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
    ggplt <- ggplot(dfplot, aes_string(x='start.genome', y='counts'))    # data
    if (both.strands) {
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='pcounts'), position=position_jitter(width=0, height=0))    # read count
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='mcounts'), position=position_jitter(width=0, height=0))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='pcounts'))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='mcounts'))    # read count
    } else {
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='counts'), position=position_jitter(width=0, height=0))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='counts'))    # read count
    }
    if (!is.null(model$segments$state)) {
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
    
    if (!is.null(model$segments$state)) {
        ggplt <- ggplt + scale_color_manual(name="state", values=stateColors(levels(dfplot.seg$state)), drop=FALSE)    # do not drop levels that are not present
    }
    if (plot.breakpoints) {
        df.bp <- as.data.frame(bp.coords)
        if (nrow(df.bp)>0) {
          statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
          suppressMessages( ggplt <- ggplt + scale_color_manual(name="state", values=c(breakpointColors(), stateColors(statelevels)), drop=FALSE) )
            ggplt <- ggplt + geom_segment(data=df.bp, aes_string(x='start.genome', xend='start.genome', color='type'), y=-1.5*custom.xlim, yend=-1.3*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'), alpha=0.5)
        }
    }
    ggplt <- ggplt + empty_theme    # no axes whatsoever
    if (both.strands) {
        ggplt <- ggplt + coord_cartesian(ylim=c(-1.5*custom.xlim,custom.xlim))    # set x- and y-limits
    } else {
        ggplt <- ggplt + coord_cartesian(ylim=c(0,custom.xlim))    # set x- and y-limits
    }
    # Get midpoints of each chromosome for xticks
    ggplt <- ggplt + scale_x_continuous(breaks=seqlengths(model$bins)/2+cum.seqlengths.0[as.character(seqlevels(model$bins))], labels=seqlevels(model$bins))
    # Quality info
    qualityInfo <- getQC(model)
    quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))
    ggplt <- ggplt + ylab(ylabstring) + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
        
    if (!is.null(file)) {
        ggsave(file, ggplt, width=20, height=5)
    } else {
        return(ggplt)
    }
}


#' Heterogeneity vs. Aneuploidy
#' 
#' Make heterogeneity vs. aneuploidy plots using individual chromosomes as datapoints.
#' 
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param hmms.list Alternative input for a faceted plot. A named list() of lists of \code{\link{aneuHMM}} objects. Alternatively a named list() of character vectors with files that contain \code{\link{aneuHMM}} objects. List names serve as facets for plotting. If specified, \code{normalChromosomeNumbers} is assumed to be a list() of vectors (or matrices) instead of a vector (or matrix).
#' @param normalChromosomeNumbers A named integer vector or matrix with physiological copy numbers, where each element (vector) or column (matrix) corresponds to a chromosome. This is useful to specify male or female samples, e.g. \code{c('X'=2)} for female samples or \code{c('X'=1,'Y'=1)} for male samples. Specify a vector if all your \code{hmms} have the same physiological copy numbers. Specify a matrix if your \code{hmms} have different physiological copy numbers (e.g. a mix of male and female samples). If not specified otherwise, '2' will be assumed for all chromosomes. If you have specified \code{hmms.list} instead of \code{hmms}, \code{normalChromosomeNumbers} is assumed to be a list() of vectors (or matrices), with one vector (or matrix) for each element in \code{hmms.list}.
#' @param plot A logical indicating whether to plot or to return the underlying data.frame.
#' @inheritParams karyotypeMeasures
#' @return A \code{\link[ggplot2]{ggplot}} object or a data.frame if \code{plot=FALSE}.
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples 
#'### Example 1: A faceted plot of lung and liver cells ###
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Make heterogeneity plots
#'plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files))
#'
#'### Example 2: Plot a mixture sample of male and female cells ###
#'## Get results from a small-cell-lung-cancer
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(lung.folder, full.names=TRUE)
#'## Construct a matrix with physiological copy numbers for a mix of 48 male and 48 female samples
#'normal.chrom.numbers <- matrix(2, nrow=96, ncol=24,
#'                               dimnames=list(sample=c(paste('male', 1:48), paste('female', 49:96)),
#'                                             chromosome=c(1:22,'X','Y')))
#'normal.chrom.numbers[1:48,c('X','Y')] <- 1
#'normal.chrom.numbers[49:96,c('Y')] <- 0
#'head(normal.chrom.numbers)
#'## Make heterogeneity plots
#'plotHeterogeneity(hmms = files, normalChromosomeNumbers = normal.chrom.numbers)
#'
#'### Example 3: A faceted plot of male lung and female liver cells ###
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Specify the physiological copy numbers
#'chrom.numbers.lung <- c(rep(2, 22), 1, 1)
#'names(chrom.numbers.lung) <- c(1:22, 'X', 'Y')
#'print(chrom.numbers.lung)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Specify the physiological copy numbers
#'chrom.numbers.liver <- c(rep(2, 22), 2, 0)
#'names(chrom.numbers.liver) <- c(1:22, 'X', 'Y')
#'print(chrom.numbers.liver)
#'## Make heterogeneity plots
#'plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files),
#'                  normalChromosomeNumbers = list(chrom.numbers.lung, chrom.numbers.liver))
#'
#'### Example 4 ###
#'## Exclude artifact regions with high variance
#'consensus <- consensusSegments(c(lung.files, liver.files))
#'variance <- apply(consensus$copy.number, 1, var)
#'exclude.regions <- consensus[variance > quantile(variance, 0.999)]
#'## Make heterogeneity plots
#'plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files),
#'                  exclude.regions=exclude.regions)
#'
plotHeterogeneity <- function(hmms, hmms.list=NULL, normalChromosomeNumbers=NULL, plot=TRUE, regions=NULL, exclude.regions=NULL) {
  
    if (!is.null(hmms.list)) {
        if (!is.null(normalChromosomeNumbers)) {
            if (!class(normalChromosomeNumbers) == 'list') {
                stop("Argument 'normalChromosomeNumbers' has to be a list with one entry (vector or matrix) for each entry in 'hmms.list'.")
            }
        }
    }
    if (is.null(hmms.list)) {
        hmms <- loadFromFiles(hmms, check.class="aneuHMM")
        ## Karyotype measures
        kmeasures <- karyotypeMeasures(hmms, normalChromosomeNumbers = normalChromosomeNumbers, regions = regions, exclude.regions = exclude.regions)
        rownames(kmeasures$genomewide) <- 'all'
        kmeasures <- rbind(kmeasures$genomewide, kmeasures$per.chromosome)
        kmeasures$chromosome <- rownames(kmeasures)
        rownames(kmeasures) <- NULL
        
        if (plot) {
            ## Plot with ggrepel
            ggplt <- ggplot(data=kmeasures, mapping=aes_string(x='Aneuploidy', y='Heterogeneity')) + geom_point()
            ggplt <- ggplt + geom_text_repel(aes_string(label='chromosome'))
            return(ggplt)
        } else {
            return(kmeasures)
        }
    } else {
        kmeasures.all <- list()
        for (i1 in 1:length(hmms.list)) {
            hmms <- hmms.list[[i1]]
            samplename <- names(hmms.list)[i1]
            hmms <- loadFromFiles(hmms, check.class="aneuHMM")
            ## Karyotype measures
            kmeasures <- karyotypeMeasures(hmms, normalChromosomeNumbers = normalChromosomeNumbers[[i1]], regions = regions, exclude.regions = exclude.regions)
            rownames(kmeasures$genomewide) <- 'all'
            kmeasures <- rbind(kmeasures$genomewide, kmeasures$per.chromosome)
            kmeasures$chromosome <- rownames(kmeasures)
            kmeasures$sample <- samplename
            kmeasures.all[[i1]] <- kmeasures
        }
        kmeasures.all <- do.call(rbind, kmeasures.all)
        rownames(kmeasures.all) <- NULL
        
        if (plot) {
            ## Plot with ggrepel
            ggplt <- ggplot(data=kmeasures.all, mapping=aes_string(x='Aneuploidy', y='Heterogeneity')) + geom_point()
            ggplt <- ggplt + geom_text_repel(aes_string(label='chromosome'))
            ggplt <- ggplt + facet_wrap(~sample)
            return(ggplt)
        } else {
            return(kmeasures.all)
        }
    }
    
}


#' Plot heatmaps for quality control
#' 
#' This function is a convenient wrapper to call \code{\link{heatmapGenomewide}} for all clusters after calling \code{\link{clusterByQuality}} and plot the heatmaps into one pdf for efficient comparison.
#' 
#' @param cl The return value of \code{\link{clusterByQuality}}.
#' @param cutree The return value of \code{\link[stats]{cutree}}, where the names correspond to the filenames to be loaded.
#' @param file A character specifying the output file.
#' @param ... Further parameters passed on to \code{\link{heatmapGenomewide}}.
#' @return A \code{\link[cowplot]{cowplot}} object or \code{NULL} if a file was specified.
#' @export
#' @examples
#'## Get a list of HMMs and cluster them
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'cl <- clusterByQuality(files, G=5)
#'heatmapGenomewideClusters(cl=cl)
#'
#'## Plot sub-clones of the largest cluster
#'largest.cluster <- which.max(sapply(cl$classification, length))
#'files <- cl$classification[[largest.cluster]]
#'clust <- clusterHMMs(files)
#'groups <- cutree(tree = clust$hclust, k = 5)
#'heatmapGenomewideClusters(cutree = groups, cluster = FALSE)
#'
heatmapGenomewideClusters <- function(cl=NULL, cutree=NULL, file=NULL, ...) {
  
    ## Check user input ##
    if (is.null(cl) & is.null(cutree)) {
        stop("Please specify either 'cl' or 'cutree'.")
    }
    if (!is.null(cl) & !is.null(cutree)) {
        stop("Please specify either 'cl' or 'cutree', not both.")
    }
  
    if (!is.null(cl)) {
        filelist <- cl$classification
    } else if (!is.null(cutree)) {
        filelist <- split(names(cutree), cutree)
    }
    ## Get the plot dimensions ##
    ptm <- startTimedMessage("Calculating plot dimensions ...")
    hmm <- loadFromFiles(filelist[[1]][1])[[1]]
      width.heatmap <- sum(as.numeric(seqlengths(hmm$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
      height <- max(length(unlist(filelist)) * 0.5, 2)
      width.dendro <- 20
      width <- width.heatmap + width.dendro
    stopTimedMessage(ptm)
  
    ## Make the plots ##
    ggplts <- list()
    for (i1 in 1:length(filelist)) {
        message("Cluster ", i1, " ...")
        ggplts[[i1]] <- heatmapGenomewide(filelist[[i1]], ...)
    }
    cowplt <- cowplot::plot_grid(plotlist = ggplts, align='v', ncol=1, rel_heights = sapply(filelist, function(x) { max(length(x), 4) }))
    if (is.null(file)) {
        return(cowplt)
    } else {
        ggsave(cowplt, filename = file, width=width, height=height, units='cm', limitsize = FALSE)
    }
    
}


#' Perform a PCA for copy number profiles
#'
#' Perform a PCA for copy number profiles in \code{\link{aneuHMM}} objects.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.
#' @param PC1 Integer specifying the first of the principal components to plot.
#' @param PC2 Integer specifying the second of the principal components to plot.
#' @param colorBy A character vector of the same length as \code{hmms} which is used to color the points in the plot.
#' @param plot Set to \code{FALSE} if you want to return the data.frame that is used for plotting instead of the plot.
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will be excluded from the computation of the PCA. This can be useful to exclude regions with artifacts.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or a data.frame if \code{plot=FALSE}.
#' @export
#' @examples
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Plot the PCA
#'classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#'labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#'plot_pca(c(lung.files, liver.files), colorBy=classes, PC1=2, PC2=4)
#' 
plot_pca <- function(hmms, PC1=1, PC2=2, colorBy=NULL, plot=TRUE, exclude.regions=NULL) {
  
    hmms <- loadFromFiles(hmms, check.class = c("aneuHMM", "aneuBiHMM"))
    copy.numbers <- sapply(hmms, function(hmm) { hmm$bins$copy.number })
    
    ### Exclude regions ###
    if (!is.null(exclude.regions)) {
        ind <- findOverlaps(hmms[[1]]$bins, exclude.regions)@from
            copy.numbers <- copy.numbers[-ind,]
    }
    
    ## PCA
    Y <- apply(log(copy.numbers+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
    s <- svd(Y)
    percent <- s$d^2/sum(s$d^2)*100
    labs <- sapply(seq_along(percent), function(i) {
        paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
    })
    df <- data.frame(PC1 = s$u[,PC1], PC2 = s$u[,PC2])
    if (!plot) {
        colnames(df) <- labs[c(PC1, PC2)]
        rownames(df) <- colnames(copy.numbers)
        return(df)
    } else {
        if (!is.null(colorBy)) {
            df$color <- colorBy
            ggplt <- ggplot(df) + geom_point(aes_string(x='PC1', y='PC2', color='color'))
            ggplt <- ggplt + scale_color_manual(values=getDistinctColors(unique(colorBy)))
        } else {
            ggplt <- ggplot(df) + geom_point(aes_string(x='PC1', y='PC2'))
        }
        ggplt <- ggplt + xlab(labs[PC1]) + ylab(labs[PC2])
        return(ggplt)
    }
}
