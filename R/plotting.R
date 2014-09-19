# ------------------------------------------------------------
# Plot a read histogram with univariate fits
# ------------------------------------------------------------
plot.distribution <- function(model, state=NULL, chrom=NULL, start=NULL, end=NULL) {

	## Load libraries
	library(ggplot2)

	# -----------------------------------------
	# Get right x limit
	get_rightxlim <- function(histdata, reads) {
		rightxlim1 <- median(reads[reads>0])*7
		breaks <- histdata$breaks[1:length(histdata$counts)]
		counts <- histdata$counts
		rightxlim2 <- breaks[counts<=5 & breaks>median(reads)*2][1]
		rightxlim <- min(c(rightxlim1,rightxlim2), na.rm=TRUE)
		return(rightxlim)
	}

	# Select the rows to plot
	selectmask <- rep(TRUE,length(model$reads))
	numchrom <- length(table(model$coordinates$chrom))
	if (!is.null(chrom)) {
		if (! chrom %in% levels(model$coordinates$chrom)) {
			stop(chrom," can't be found in the model coordinates.")
		}
		selectchrom <- model$coordinates$chrom == chrom
		selectmask <- selectmask & selectchrom
		numchrom <- 1
	}
	if (numchrom == 1) {
		if (!is.null(start)) {
			selectstart <- model$coordinates$start >= start
			selectmask <- selectmask & selectstart
		}
		if (!is.null(end)) {
			selectend <- model$coordinates$end <= end
			selectmask <- selectmask & selectend
		}
	}
	if (!is.null(state)) {
		selectmask <- selectmask & model$states==state
	}
	if (length(which(selectmask)) != length(model$reads)) {
		reads <- model$reads[selectmask]
# 		posteriors <- model$posteriors[selectmask,]
# 		weights <- apply(posteriors,2,mean)
		states <- model$states[selectmask]
		weights <- rep(NA, 3)
		for (i1 in model$use.states+1) {
			weights[i1] <- length(which(states==state.labels[i1]))
		}
		weights <- weights / length(states)
	} else {
		reads <- model$reads
		weights <- model$weights
	}

	# Find the x limits
	breaks <- max(reads)
	if (max(reads)==0) { breaks <- 1 }
	histdata <- hist(reads, right=FALSE, breaks=breaks, plot=FALSE)
	rightxlim <- get_rightxlim(histdata, reads)

	# Plot the histogram
	ggplt <- ggplot(data.frame(reads)) + geom_histogram(aes(x=reads, y=..density..), binwidth=1, color='black', fill='white') + xlim(0,rightxlim) + theme_bw() + xlab("read count")

	### Add fits to the histogram
	numstates <- model$num.states
	x <- 0:rightxlim
	distributions <- list(x)

	# zero-inflation
	distributions[[length(distributions)+1]] <- c(weights[1],rep(0,length(x)-1))
	# geometric
	distributions[[length(distributions)+1]] <- weights[2] * dgeom(x, model$distributions[2,'prob'])
	# negative binomials
	for (istate in 3:numstates) {
		distributions[[length(distributions)+1]] <- weights[istate] * dnbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])
	}
	distributions <- as.data.frame(distributions)
	names(distributions) <- c("x",state.labels[model$use.states+1])
	# Total
	distributions$total <- apply(distributions[-1], 1, sum)

	# Reshape the data.frame for plotting with ggplot
	distributions <- reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="state", times=c(state.labels[model$use.states+1],"total"))
	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=state, col=state), data=distributions)
	} else {
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=state, size=state), data=distributions[c(1,state+1)])
	}
	
	# Make legend and colors correct
	lmeans <- round(model$distributions[,'mu'], 2)
	lvars <- round(model$distributions[,'variance'], 2)
	legend <- paste(state.labels, ", mean=", lmeans, ", var=", lvars, sep='')
	legend <- c(legend, paste0('total, mean(data)=', round(mean(reads),2), ', var(data)=', round(var(reads),2)))
	ggplt <- ggplt + scale_color_manual(breaks=c(state.labels, 'total'), values=state.colors[c(state.labels,'total')], labels=legend)
	ggplt <- ggplt + theme(legend.position=c(1,1), legend.justification=c(1,1))

	return(ggplt)

}


# ------------------------------------------------------------
# Plot a boxplot of the univariate calls
# ------------------------------------------------------------
plot.boxplot <- function(model) {

	## Load libraries
	library(ggplot2)

	## Plot settings
	cols <- gcolors[c("unmodified","modified")]

	## Boxplot
	components <- c("unmodified","modified")[as.factor(get.states(model))]
	df <- data.frame(component=components, reads=model$reads)
	ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes(x=component, y=reads, fill=component)) + scale_fill_manual(values=cols)
	return(ggplt)

}


# ------------------------------------------------------------
# Plot like BAIT
# ------------------------------------------------------------
plot.BAIT <- function(model, file='aneufinder_BAIT_plots') {
	
	## Convert to GRanges
	gr <- hmm2GRanges(model, reduce=F)
	grl <- split(gr, seqnames(gr))

	## Get some variables
	num.chroms <- length(levels(seqnames(gr)))
	maxseqlength <- max(seqlengths(gr))
	custom.xlim <- model$distributions['monosomy','mu'] * 10

	## Setup page
	library(grid)
	library(ggplot2)
	nrows <- 2	# rows for plotting chromosomes
	ncols <- ceiling(num.chroms/nrows)
	png(file=paste0(file, '.png'), width=ncols*6, height=nrows*21, units='cm', res=150)
	grid.newpage()
	layout <- matrix(1:((nrows+1)*ncols), ncol=ncols, nrow=nrows+1, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,21,21))))
	# Main title
	grid.text(file, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=26))

	## Go through chromosomes and plot
	for (i1 in 1:num.chroms) {
		# Get the i,j matrix positions of the regions that contain this subplot
		matchidx <- as.data.frame(which(layout == i1+ncols, arr.ind = TRUE))

		# Percentage of chromosome in state
		tstate <- table(mcols(grl[[i1]])$state)
		pstate.all <- tstate / sum(tstate)
		pstate <- round(pstate.all*100)[-1]	# without 'nullsomy / unmapped' state
		pstring <- apply(pstate, 1, function(x) { paste0(": ", x, "%") })
		pstring <- paste0(names(pstring), pstring)
		pstring <- paste(pstring[which.max(pstate)], collapse="\n")
		pstring2 <- round(pstate.all*100)[1]	# only 'nullsomy / unmapped'
		pstring2 <- paste0(names(pstring2), ": ", pstring2, "%")

		# Plot the read counts
		dfplot <- as.data.frame(grl[[i1]])
		dfplot.points <- dfplot[dfplot$reads>=custom.xlim,]
		dfplot$reads <- stats::runmed(dfplot$reads, 11)
		empty_theme <- theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
			axis.title.x=element_text(size=20),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
		ggplt <- ggplot(dfplot, aes(x=start, y=reads))	# data
		ggplt <- ggplt + geom_linerange(aes(ymin=0, ymax=reads, col=state), size=0.2)	# read count
		ggplt <- ggplt + xlim(0,maxseqlength) + ylim(-0.6*custom.xlim,custom.xlim)	# set x- and y-limits
		ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, mapping=aes(xmax=max(start)), col='white', fill='gray20')	# chromosome backbone as simple rectangle
		ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes(x=start, y=reads, col=state))	# outliers
		ggplt <- ggplt + scale_color_manual(values=state.colors, drop=F)	# do not drop levels if not present
		ggplt <- ggplt + empty_theme	# no axes whatsoever
		ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1], "\n", pstring, "\n", pstring2))	# chromosome names
		ggplt <- ggplt + coord_flip()
		suppressWarnings(
		print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		)
		
	}
	d <- dev.off()
}



# ------------------------------------------------------------
# Plot overview
# ------------------------------------------------------------
plot.genome.overview <- function(modellist, file='aneufinder_genome_overview', numCPU=1) {
	
	## Function definitions
	reformat <- function(x) {
	out_list <- list() 

		for ( i in 2:length(x) ) {
			out_list[[i]] <- c(x[i-1], x[i])
		}
	mt <- do.call("rbind",out_list)
	df <- data.frame(mt)
	colnames(df) <- c("start", "end")
	df
	}

	## Load and transform to GRanges
	cat('transforming to GRanges\n')
	uni.hmm.grl <- hmmList2GRangesList(modellist, reduce=FALSE)

	## Setup page
	library(grid)
	library(ggplot2)
	nrows <- length(uni.hmm.grl)	# rows for plotting genomes
	ncols <- 1
	png(file=paste0(file, '.png'), width=ncols*60, height=nrows*5, units='cm', res=150)
	grid.newpage()
	layout <- matrix(1:((nrows+1)*ncols), ncol=ncols, nrow=nrows+1, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,rep(10,length(uni.hmm.grl))))))
	# Main title
	grid.text(file, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=26))

	## Prepare some variables for plotting
	gr <- uni.hmm.grl[[1]]
	len <- seqlengths(gr)
	chr.names <- levels(seqnames(gr))
	len <- as.numeric(len)
	len <- c(0,len)
	len <- cumsum(len)
	df <- reformat(len)
	df$col <- rep(c("grey47","grey77"), 12)
	df$breaks <- df[,1] + ((df[,2]-df[,1])/2)
	df$chr.names <- chr.names 

	my_theme <- theme(
				legend.position="none",
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
	)
		
	## Go through models and plot
	for (i1 in 1:length(uni.hmm.grl)) {
		cat('plotting model',i1,'\n')
		# Get the i,j matrix positions of the regions that contain this subplot
		matchidx <- as.data.frame(which(layout == i1+ncols, arr.ind = TRUE))

		trans_gr <- biovizBase::transformToGenome(uni.hmm.grl[[i1]], space.skip = 0)

		dfplot <- as.data.frame(trans_gr)
		dfplot$reads <- stats::runmed(uni.hmm.grl[[i1]]$reads, 15)

		ggplt <- ggplot(dfplot, aes(x=.start, y=reads))
		ggplt <- ggplt + geom_linerange(aes(ymin=0, ymax=reads, col=state), size=0.2)
		ggplt <- ggplt + geom_rect(data=df, aes(xmin=start, xmax=end, ymin=0, ymax=Inf, fill = col),  alpha=I(0.3), inherit.aes = F)
		ggplt <- ggplt + scale_x_continuous(breaks = df$breaks, labels=df$chr.names, expand = c(0,0)) + scale_y_continuous(expand=c(0,5)) + scale_fill_manual(values = c("grey47","grey77")) + scale_color_manual(values=state.colors, drop=F) + xlab("chromosomes") + my_theme
		suppressWarnings(
		print(ggplt + ylim(0,30), vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		)


	}
	d <- dev.off()
}

