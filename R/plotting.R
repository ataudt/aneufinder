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
		rightxlim2 <- breaks[counts<=5 & breaks>median(reads)][1]
		rightxlim <- min(c(rightxlim1,rightxlim2), na.rm=TRUE)
		return(rightxlim)
	}

	## Plot settings
	cols <- gcolors[c("unmappaple","monosomy","disomy","trisomy","tetrasomy")]

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

	# unmappable
	distributions[[length(distributions)+1]] <- c(weights[1],rep(0,length(x)-1))
	for (istate in 2:numstates) {
		distributions[[length(distributions)+1]] <- weights[istate] * dnbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])
	}
	distributions <- as.data.frame(distributions)
	names(distributions) <- c("x",state.labels[model$use.states+1])
	# Total
	distributions$total <- apply(distributions[-1], 1, sum)

	# Reshape the data.frame for plotting with ggplot
	distributions <- reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="xsomy", times=c(state.labels[model$use.states+1],"total"))
	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=xsomy, cols=xsomy), data=distributions[distributions$xsomy!="total",])
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=xsomy), data=distributions[distributions$xsomy=="total",])
	} else {
		ggplt <- ggplt + geom_line(aes(x=x, y=density, group=xsomy, size=xsomy), data=distributions[c(1,state+1)])
	}
	
	# Make legend and colors correct
	ggplt <- ggplt + scale_color_manual(name="components", values=cols)

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

