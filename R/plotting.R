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
		for (i1 in 1:length(levels(model$states))) {
			weights[i1] <- length(which(states==model$state.labels[i1]))
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
	c.state.labels <- as.character(model$state.labels)
	numstates <- model$num.states
	x <- 0:rightxlim
	distributions <- list(x)

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
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=state, col=state), data=distributions)
	} else {
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=state, size=state), data=distributions[distributions$state==state,])
	}
	
	# Make legend and colors correct
	lmeans <- round(model$distributions[,'mu'], 2)
	lvars <- round(model$distributions[,'variance'], 2)
	legend <- paste(c.state.labels, ", mean=", lmeans, ", var=", lvars, sep='')
	legend <- c(legend, paste0('total, mean(data)=', round(mean(reads),2), ', var(data)=', round(var(reads),2)))
	ggplt <- ggplt + scale_color_manual(breaks=c(c.state.labels, 'total'), values=state.colors[c(c.state.labels,'total')], labels=legend)
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
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.chromosomes <- function(model, file=NULL) {

	if (class(model)==class.univariate.hmm) {
		plot.chromosomes.univariate(model, file=file)
	} else if (class(model)==class.multivariate.hmm) {
		plot.chromosomes.bivariate(model, file=file)
	}

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.chromosomes.univariate <- function(model, file=NULL) {
	
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
	if (!is.null(file)) {
		pdf(file=paste0(file, '.pdf'), width=ncols*2.4, height=nrows*8.3)
	}
	grid.newpage()
	layout <- matrix(1:((nrows+1)*ncols), ncol=ncols, nrow=nrows+1, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,21,21))))
	# Main title
	grid.text(model$ID, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=26))

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
		# Set values to great for plotting to limit
			dfplot$reads[dfplot$reads>=custom.xlim] <- custom.xlim
			dfplot.points <- dfplot[dfplot$reads>=custom.xlim,]
			dfplot.points$reads <- rep(custom.xlim, nrow(dfplot.points))
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
		ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, mapping=aes(xmax=max(start)), col='white', fill='gray20')	# chromosome backbone as simple rectangle
		ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes(x=start, y=reads, col=state), size=5, shape=21)	# outliers
		ggplt <- ggplt + scale_color_manual(values=state.colors, drop=F)	# do not drop levels if not present
		ggplt <- ggplt + empty_theme	# no axes whatsoever
		ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1], "\n", pstring, "\n", pstring2))	# chromosome names
		ggplt <- ggplt + xlim(0,maxseqlength) + ylim(-0.6*custom.xlim,custom.xlim)	# set x- and y-limits
		ggplt <- ggplt + coord_flip()
		suppressWarnings(
		print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		)
		
	}
	if (!is.null(file)) {
		d <- dev.off()
	}
}



# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.chromosomes.bivariate <- function(model, file=NULL) {
	
	## Convert to GRanges
	gr <- hmm2GRanges(model, reduce=F)
	grl <- split(gr, seqnames(gr))

	## Get some variables
	num.chroms <- length(levels(seqnames(gr)))
	maxseqlength <- max(seqlengths(gr))
	custom.xlim <- model$distributions.univariate[[1]]['monosomy','mu'] * 10

	## Setup page
	library(grid)
	library(ggplot2)
	nrows <- 2	# rows for plotting chromosomes
	ncols <- ceiling(num.chroms/nrows)
	if (!is.null(file)) {
		pdf(file=paste0(file, '.pdf'), width=ncols*2.4, height=nrows*8.3)
	}
	grid.newpage()
	layout <- matrix(1:((nrows+1)*ncols), ncol=ncols, nrow=nrows+1, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,21,21))))
	# Main title
	grid.text(model$ID, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=26))

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

		## Convert to data.frame for plotting and prepare the data.frame
			dfplot <- as.data.frame(grl[[i1]])
			dfplot$reads.minus <- - dfplot$reads.minus
		# Set values to great for plotting to limit
			dfplot$reads.plus[dfplot$reads.plus>=custom.xlim] <- custom.xlim
			dfplot$reads.minus[dfplot$reads.minus<=-custom.xlim] <- -custom.xlim
			dfplot.points.plus <- dfplot[dfplot$reads.plus>=custom.xlim,]
			dfplot.points.plus$reads <- rep(custom.xlim, nrow(dfplot.points.plus))
			dfplot.points.minus <- dfplot[dfplot$reads.minus<=-custom.xlim,]
			dfplot.points.minus$reads <- rep(-custom.xlim, nrow(dfplot.points.minus))
		# No axes and labels
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
		# Plot
		ggplt <- ggplot(dfplot, aes(x=start, y=reads.plus))	# data
		ggplt <- ggplt + geom_linerange(aes(ymin=0, ymax=reads.plus, col=state.separate.plus), size=0.2)	# read count
		ggplt <- ggplt + geom_linerange(aes(ymin=0, ymax=reads.minus, col=state.separate.minus), size=0.2)	# read count
		ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=0, mapping=aes(xmax=max(start)), col='white', fill='gray20')	# chromosome backbone as simple rectangle
		ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes(x=start, y=reads, col=state.separate.plus), size=5, shape=21)	# outliers
		ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes(x=start, y=reads, col=state.separate.minus), size=5, shape=21)	# outliers
		ggplt <- ggplt + scale_color_manual(values=state.colors, drop=F)	# do not drop levels if not present
		ggplt <- ggplt + empty_theme	# no axes whatsoever
		ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1], "\n", pstring, "\n", pstring2))	# chromosome names
		ggplt <- ggplt + xlim(0,maxseqlength) + ylim(-custom.xlim,custom.xlim)	# set x- and y-limits
		ggplt <- ggplt + coord_flip()
		suppressWarnings(
		print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		)
		
	}
	if (!is.null(file)) {
		d <- dev.off()
	}
}



# ------------------------------------------------------------
# Plot overview
# ------------------------------------------------------------
plot.genome.overview <- function(hmm.list, file='aneufinder_genome_overview', numCPU=1, chromosome='') {
	
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

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)
	

	## Load and transform to GRanges
	cat('transforming to GRanges\n')
	temp <- hmmList2GRangesList(hmm.list, reduce=FALSE, numCPU=numCPU)
	uni.hmm.grl <- temp$grl

	## Setup page
	library(grid)
	library(ggplot2)
	nrows <- length(uni.hmm.grl)	# rows for plotting genomes
	ncols <- 1
	pdf(file=paste0(file, '.pdf'), width=ncols*24, height=nrows*2)
	grid.newpage()
	layout <- matrix(1:((nrows+1)*ncols), ncol=ncols, nrow=nrows+1, byrow=T)
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,rep(10,length(uni.hmm.grl))))))
	# Main title
	grid.text(file, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=26))

	## Prepare some variables for plotting
	gr <- uni.hmm.grl[[1]]
	len <- seqlengths(gr)
	chr.names <- names(seqlengths(gr))
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

	##get the max read count in GRangesList
	max_list <- lapply(uni.hmm.grl, function(x) max(stats::runmed(x$reads, 15)))  #stats::runmed = filtering extreme read counts values
	max_reads <- max(unlist(max_list))
		
	## Go through models and plot
	for (i1 in 1:length(uni.hmm.grl)) {
		cat('plotting model',i1,'\n')
		# Get the i,j matrix positions of the regions that contain this subplot
		matchidx <- as.data.frame(which(layout == i1+ncols, arr.ind = TRUE))

		trans_gr <- biovizBase::transformToGenome(uni.hmm.grl[[i1]], space.skip = 0)

		dfplot <- as.data.frame(trans_gr)

		ggplt <- ggplot(dfplot, aes(x=.start, y=reads))
		ggplt <- ggplt + geom_linerange(aes(ymin=0, ymax=reads, col=state), size=0.2)
		ggplt <- ggplt + geom_rect(data=df, aes(xmin=start, xmax=end, ymin=0, ymax=Inf, fill = col),  alpha=I(0.3), inherit.aes = F)
		if (!chromosome == '') { ##zoom into the specific chromosome
			xlim <- df[chromosome,]
			ggplt <- ggplt + scale_x_continuous(breaks = df$breaks, labels=df$chr.names, expand = c(0,0)) + scale_y_continuous(expand=c(0,5)) + scale_fill_manual(values = c("grey47","grey77")) + scale_color_manual(values=state.colors, drop=F) + xlab("chromosome") + my_theme + coord_cartesian(xlim = c(xlim$start, xlim$end))
			suppressWarnings(
			print(ggplt + ylim(0,max_reads), vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
			)
		} else {
			ggplt <- ggplt + scale_x_continuous(breaks = df$breaks, labels=df$chr.names, expand = c(0,0)) + scale_y_continuous(expand=c(0,5)) + scale_fill_manual(values = c("grey47","grey77")) + scale_color_manual(values=state.colors, drop=F) + xlab("chromosomes") + my_theme
			suppressWarnings(
			print(ggplt + ylim(0,max_reads), vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
			)
		}
	}
	d <- dev.off()
}


# ------------------------------------------------------------
# Plot a heatmap of chromosome state for multiple samples
# ------------------------------------------------------------
plot.chromosome.heatmap <- function(hmm.list, cluster=TRUE, numCPU=1) {

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)
	
	## Transform to GRanges in reduced representation
	temp <- hmmList2GRangesList(hmm.list, reduce=T, numCPU=numCPU, consensus=F)
	grlred <- temp$grl
	
	## Find the most frequent state (mfs) for each chromosome and sample
	cat("finding most frequent state for each sample and chromosome ...")
	grl.per.chrom <- lapply(grlred, function(x) { split(x, seqnames(x)) })
	mfs.samples <- list()
	for (i1 in 1:length(hmm.list)) {
		mfs.samples[[hmm.list[[i1]]$ID]] <- lapply(grl.per.chrom[[i1]], function(x) { tab <- aggregate(width(x), by=list(state=x$state), FUN="sum"); tab$state[which.max(tab$x)] })
		attr(mfs.samples[[hmm.list[[i1]]$ID]], "varname") <- 'chromosome'
	}
	attr(mfs.samples, "varname") <- 'sample'
	cat(" done\n")

	## Transform to data.frame
	# Long format
	df <- reshape2::melt(mfs.samples, value.name='state')
	df$state <- factor(df$state, levels=levels(hmm.list[[1]]$states))
	df$sample <- factor(df$sample, levels=unique(df$sample))
	df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))

	## Cluster the samples by chromosome state
	if (cluster) {
		# Wide format
		df.wide <- reshape2::dcast(df, sample ~ chromosome, value.var='state', factorsAsStrings=F)
		# Correct strings to factors
		for (col in 2:ncol(df.wide)) {
			df.wide[,col] <- factor(df.wide[,col], levels=levels(hmm.list[[1]]$states))
		}
		# Cluster
		hc <- hclust(dist(data.matrix(df.wide[-1])))
		# Reorder samples in mfs list
		mfs.samples.clustered <- mfs.samples[hc$order]
		attr(mfs.samples.clustered, "varname") <- 'sample'
		df <- reshape2::melt(mfs.samples.clustered, value.name='state')
		df$state <- factor(df$state, levels=levels(hmm.list[[1]]$states))
		df$sample <- factor(df$sample, levels=unique(df$sample))
		df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
	}

	## Plot to heatmap
	ggplt <- ggplot(df) + geom_tile(aes(x=chromosome, y=sample, fill=state), col='black') + theme_bw() + scale_fill_manual(values=get.state.colors()[levels(df$state)])
	return(ggplt)
}


# ------------------------------------------------------------
# Plot a clustered heatmap of state calls
# ------------------------------------------------------------
plot.genome.heatmap <- function(hmm.list, file=NULL, cluster=TRUE, numCPU=1) {

	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)

	if (cluster) {
		# Transform to GRanges in reduced representation
		temp <- hmmList2GRangesList(hmm.list, reduce=TRUE, numCPU=numCPU, consensus=T)
		grlred <- temp$grl
		consensus <- temp$consensus
		constates <- temp$constates
		# Distance measure
		cat("clustering ...")
		constates[is.na(constates)] <- 0
		wcor <- cov.wt(constates, wt=as.numeric(width(consensus)), cor=T)
		dist <- as.dist(1-wcor$cor)
		# Dendrogram
		hc <- hclust(dist)
		# Reorder samples
		grlred <- grlred[hc$order]
		cat(" done\n")
	} else {
		# Transform to GRanges in reduced representation
		temp <- hmmList2GRangesList(hmm.list, reduce=TRUE, numCPU=numCPU, consensus=F)
		grlred <- temp$grl
	}

	cat("transforming coordinates ...")
	## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
	cum.seqlengths <- cumsum(as.numeric(seqlengths(grlred[[1]])))
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- seqlevels(grlred[[1]])
	transCoord <- function(gr) {
		gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
		gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
		return(gr)
	}
	grlred <- endoapply(grlred, transCoord)
	cat(" done\n")

	## Data.frame for plotting
	cat("making the plot ...")
	# Data
	df <- list()
	for (i1 in 1:length(grlred)) {
		df[[length(df)+1]] <- data.frame(start=grlred[[i1]]$start.genome, end=grlred[[i1]]$end.genome, seqnames=seqnames(grlred[[i1]]), sample=hmm.list[[i1]]$ID, state=grlred[[i1]]$state)
	}
	df <- do.call(rbind, df)
	# Chromosome lines
	label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grlred[[1]]) )
	df.chroms <- data.frame(y=c(0,cum.seqlengths))

	## Plot
	ggplt <- ggplot(df) + geom_linerange(aes(ymin=start, ymax=end, x=sample, col=state), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=state.colors) + theme(panel.background=element_blank(), axis.ticks.x=element_blank())
	ggplt <- ggplt + geom_hline(aes(yintercept=y), data=df.chroms, col='black')
	cat(" done\n")

	## Plot to file
	if (!is.null(file)) {
		cat(paste0("plotting to file ",file," ..."))
		height.cm <- length(hmm.list) * 0.5
		width.cm <- 200
		pdf(file, width=width.cm/2.54, height=height.cm/2.54)
		print(ggplt)
		d <- dev.off()
		cat(" done\n")
	} else {
		return(ggplt)
	}

}
