# ------------------------------------------------------------
# Plot a read histogram with univariate fits
# ------------------------------------------------------------
plot.distribution <- function(model, state=NULL, chr=NULL, start=NULL, end=NULL) {

	## Load libraries
# 	library(ggplot2)

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
	cols <- c("unmodified"="gray48","modified"="orangered3", "total"="black")

	# Select the rows to plot
	selectmask <- rep(TRUE,length(model$reads))
	numchrom <- length(table(model$coordinates$chrom))
	if (!is.null(chr)) {
		if (! chr %in% levels(model$coordinates$chrom)) {
			stop(chr," can't be found in the model coordinates.")
		}
		selectchrom <- model$coordinates$chrom == chr
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
		states <- get.states(model)
		selectmask <- selectmask & states==state
	}
	if (length(which(selectmask)) != length(model$reads)) {
		reads <- model$reads[selectmask]
		posteriors <- model$posteriors[selectmask,]
		softweights <- apply(posteriors,2,mean)
	} else {
		reads <- model$reads
		softweights <- model$softweights
	}

	# Find the x limits
	breaks <- max(reads)
	if (max(reads)==0) { breaks <- 1 }
	histdata <- hist(reads, right=FALSE, breaks=breaks, plot=FALSE)
	rightxlim <- get_rightxlim(histdata, reads)

	# Plot the histogram
	ggplt <- ggplot(data.frame(reads)) + geom_histogram(aes(x=reads, y=..density..), binwidth=1, color='black', fill='white') + xlim(0,rightxlim) + theme_bw() + xlab("read count")

	### Add fits to the histogram
	numstates <- length(softweights)
	x <- 0:rightxlim
	distributions <- data.frame(x)

	# Unmodified
	distributions$unmodified <- (1-softweights[3]) * dzinbinom(x, softweights[1], model$distributions[2,'r'], model$distributions[2,'p'])
	# Modified
	distributions$modified <- softweights[3] * dnbinom(x, model$distributions[3,'r'], model$distributions[3,'p'])
	# Total
	distributions$total <- distributions$unmodified + distributions$modified

	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(aes(x=x, y=unmodified, color="unmodified", group=1), data=distributions, size=1)
		ggplt <- ggplt + geom_line(aes(x=x, y=modified, color="modified", group=1), data=distributions, size=1)
		ggplt <- ggplt + geom_line(aes(x=x, y=total, color="total", group=1), data=distributions, size=1)
	} else {
		if (state==0) ggplt <- ggplt + geom_line(aes(x=x, y=unmodified, color="unmodified", group=1), data=distributions, size=1)
		if (state==1) ggplt <- ggplt + geom_line(aes(x=x, y=modified, color="modified", group=1), data=distributions, size=1)
	}
	
	# Make legend and colors correct
	ggplt <- ggplt + scale_color_manual(name="components", values=cols)

	return(ggplt)

}

# ------------------------------------------------------------
# Plot a read histogram in normal space of the given state
# ------------------------------------------------------------
plot.distribution.normal <- function(model, state=0) {

	## Load libraries
# 	library(ggplot2)

	## Intercept user input
	if (state!=0 & state!=1) { stop("state has to be either 0 or 1") }

	## Plot settings
	cols <- c("unmodified"="gray48","modified"="orangered3")

	## Transform the reads
	states <- get.states(model)
	df <- data.frame(bin=1:length(model$reads), reads=model$reads, state=as.factor(states))
	# Transform to uniform space
	df$ureads[df$state==0] <- pzinbinom(df$reads[df$state==0], model$softweights[1], model$distributions[2,'r'], model$distributions[2,'p'])
	df$ureads[df$state==1] <- pnbinom(df$reads[df$state==1], model$distributions[3,'r'], model$distributions[3,'p'])
	# Transform to normal space
	df$nreads <- qnorm(df$ureads)

	## Make the plots
	subset <- df$nreads[df$state==state]
	breaks <- c(-Inf,sort(as.numeric(names(table(subset)))))
	x <- seq(-4,4,0.1)
	title <- paste("Transformed emission density for state ",state, sep="")
	ggplt <- ggplot() + geom_histogram(data=data.frame(ureads=subset), aes(x=ureads, y=..density..), breaks=breaks, right=TRUE, col='black', fill=cols[state+1]) + theme_bw() + geom_line(data=data.frame(x=x, y=dnorm(x, mean=0, sd=1)), aes(x=x, y=y)) + xlab("transformed reads") + ylim(0,0.5) + labs(title=title)
	return(ggplt)

}

# ------------------------------------------------------------
# Plot a boxplot of the univariate calls
# ------------------------------------------------------------
plot.boxplot <- function(model) {

	## Load libraries
# 	library(ggplot2)

	## Plot settings
	cols <- c("unmodified"="gray48","modified"="orangered3")

	## Boxplot
	components <- c("unmodified","modified")[as.factor(get.states(model))]
	df <- data.frame(component=components, reads=model$reads)
	ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes(x=component, y=reads, fill=component)) + scale_fill_manual(values=cols)
	return(ggplt)

}

