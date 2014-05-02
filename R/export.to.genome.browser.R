# ----------------------------------------------------------------
# Write univariate states to file
# ----------------------------------------------------------------
univariate2bed <- function(modellist, file="view_me_in_genome_browser", threshold=0.5, separate.zeroinflation=FALSE, chrom.length.file=NULL) {

	# Check user input
	if (is.null(names(modellist))) {
		stop("Please name the list entries. The names will be used as track name for the Genome Browser file.")
	}

	# Variables
	nummod <- length(modellist)
	file <- paste(file,"bed", sep=".")

	# Generate the colors
	colors <- c("zeroinflation"="black", "unmodified"="gray48","modified"="orangered3")
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")

	# Write first line to file
	cat("browser hide all\n", file=file)
	
	for (imod in 1:nummod) {
		model <- modellist[[imod]]
		priority <- 50 + imod
		cat(paste("track name=",names(modellist)[imod]," description=\"univariate calls for ",names(modellist)[imod],", threshold = ",threshold,"\" visibility=1 itemRgb=On priority=",priority,"\n", sep=""), file=file, append=TRUE)

		# Collapse the calls
		if (separate.zeroinflation) {
			calls <- cbind(model$coordinates, name=model$states.with.zeroinflation)
			collapsed.calls <- collapse.bins(calls, column2collapseBy=4)
		} else {
			calls <- cbind(model$coordinates, name=model$states)
			collapsed.calls <- collapse.bins(calls, column2collapseBy=4)
			collapsed.calls <- collapsed.calls[collapsed.calls$name=="modified",]
		}

		# Check length of chromosomes if chrom.length.file was given
		if (!is.null(chrom.length.file)) {
			chrom.lengths.df <- read.table(chrom.length.file)
			chrom.lengths <- chrom.lengths.df[,2]
			names(chrom.lengths) <- chrom.lengths.df[,1]
			for (chrom in levels(as.factor(collapsed_calls$chrom))) {
				index <- which(collapsed_calls$chrom == chrom & collapsed_calls$end > chrom.lengths[chrom])
				collapsed_calls$end[index] <- chrom.lengths[chrom]
				if (length(index) >= 1) {
					warning("Adjusted entry in chromosome ",chrom,", because it was longer than the length of the chromosome")
				}
			}
		}

		# Append to file
		itemRgb <- RGBs[as.character(collapsed.calls$name)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		write.table(format(df, scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

}

# ----------------------------------------------------------------
# 
# ----------------------------------------------------------------
make.ucsc.file.combinatorial <- function(modellist, numstates=NULL, file="combinatorial-states-ucsc", threshold=0.5, chrom.length.file=NULL) {
	
	# Define variables
	nummod <- length(modellist)
	file <- paste(file,".BED", sep="")
	coordinates <- modellist[[1]]$coordinates
	
	# Write first line to file
	cat("browser position chr1:1-100000\n", file=file)

	# Get the combinatorial states
	cat("Getting the combinatorial states ...\n")
	combstates <- combinatorial.states(modellist)
	states <- as.numeric(names(sort(table(combstates), decreasing=TRUE)))

	# Clean up to reduce memory usage
	remove(modellist)

	if (is.null(numstates)) {
		numstates <- length(states)
	} else if (numstates > length(states)) {
		numstates <- length(states)
	}
	# Go through each state and write to file
	for (istate in states[1:numstates]) {
		cat("\nAnalyzing state",istate,"\n")
		cat(paste("track name=\"combinatorial state ",istate,"\" description=\"univariate calls for state ",istate,", threshold = ",threshold,"\" visibility=1\n", sep=""), file=file, append=TRUE)

		calls <- coordinates[ combstates == istate , ]
		if (length(calls$start) >= 1) {
			# Collapse the bins
			cat("collapsing the calls\n")
			collapsed_calls <- collapse.bins(calls)
		} else {
			cat("no bins in state",istate,"\n")
			collapsed_calls <- NULL
		}

		if (!is.null(collapsed_calls)) {
			# Check length of chromosomes if chrom.length.file was given
			if (!is.null(chrom.length.file)) {
				chrom.lengths.df <- read.table(chrom.length.file)
				chrom.lengths <- chrom.lengths.df[,2]
				names(chrom.lengths) <- chrom.lengths.df[,1]
				for (chrom in levels(as.factor(collapsed_calls$chrom))) {
					index <- which(collapsed_calls$chrom == chrom & collapsed_calls$end > chrom.lengths[chrom])
					collapsed_calls$end[index] <- chrom.lengths[chrom]
					if (length(index) >= 1) {
						warning("Adjusted entry in chromosome ",chrom,", because it was longer than the length of the chromosome")
					}
				}
			}

			cat("appending to file ...\n")
			write.table(format(as.data.frame(collapsed_calls)[,1:3], scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}

}


# ----------------------------------------------------------------
# 
# ----------------------------------------------------------------
make.ucsc.file.multivariate <- function(model, numstates=NULL, file="multivariate-states-ucsc", chrom.length.file=NULL) {

	# Variable assignments
	states <- model$stateorder
	coordinates <- model$coordinates
	file <- paste(file,".BED", sep="")

	# Write first line to file
	cat("browser position chr1:1-100000\n", file=file)

	# Determine state of bin
	cat("Determining state of bins ...\n")
	combstates <- apply(model$posteriors,1,which.max)
	combstates <- states[combstates]

	if (is.null(numstates)) {
		numstates <- length(states)
	} else if (numstates > length(states)) {
		numstates <- length(states)
	}
	# Go through each state and write to file
	for (istate in states[1:numstates]) {
		cat("\nAnalyzing state",istate,"\n")
		cat(paste("track name=\"multivariate state ",istate,"\" description=\"multivariate calls for state ",istate,"\" visibility=1 color=0,0,255\n", sep=""), file=file, append=TRUE)

		calls <- coordinates[ combstates == istate , ]
		if (length(calls$start) >= 1) {
			# Collapse the bins
			cat("collapsing the calls\n")
			collapsed_calls <- collapse.bins(calls)
		} else {
			cat("no bins in state",istate,"\n")
			collapsed_calls <- NULL
		}

		if (!is.null(collapsed_calls)) {
			# Check length of chromosomes if chrom.length.file was given
			if (!is.null(chrom.length.file)) {
				chrom.lengths.df <- read.table(chrom.length.file)
				chrom.lengths <- chrom.lengths.df[,2]
				names(chrom.lengths) <- chrom.lengths.df[,1]
				for (chrom in levels(as.factor(collapsed_calls$chrom))) {
					index <- which(collapsed_calls$chrom == chrom & collapsed_calls$end > chrom.lengths[chrom])
					collapsed_calls$end[index] <- chrom.lengths[chrom]
					if (length(index) >= 1) {
						warning("Adjusted entry in chromosome ",chrom,", because it was longer than the length of the chromosome")
					}
				}
			}

			cat("appending to file ...\n")
			write.table(format(as.data.frame(collapsed_calls)[,1:3], scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}

}


# ----------------------------------------------------------------
# Write signal tracks from binned data
# ----------------------------------------------------------------
bin2wiggle <- function(binned.data.list, file="view_me_in_genome_browser") {
	
	# Check user input
	if (is.null(names(binned.data.list))) {
		stop("Please name the list entries. The names will be used as track name for the Genome Browser file.")
	}

	# Variables
	file <- paste(file,"wig", sep=".")

	# Write to file
	cat("browser hide all\n", file=file)
	
	# Go through binned.data.list
	for (i1 in 1:length(binned.data.list)) {
		binned.data <- binned.data.list[[i1]]
		names(binned.data) <- c("chrom","start","end","reads")
		
		# Write track information
		name <- names(binned.data.list)[i1]
		binsize <- diff(binned.data$start[1:2])
		description <- paste("binned read counts ",binsize,"bp", sep="")
		cat(paste("track type=wiggle_0 name=\"",name,"\" description=\"",description,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",i1,"\n", sep=""), file=file, append=TRUE)
		# Write read data
		for (chrom in unique(binned.data$chrom)) {
			cat(paste("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n", sep=""), file=file, append=TRUE)
			write.table(binned.data$reads, file=file, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}

}


# ----------------------------------------------------------------
# Write color-coded (multivariate) combinatorial states to file
# ----------------------------------------------------------------
multivariate2bed <- function(model, separate.tracks=FALSE, exclude.state.zero=TRUE, numstates=NULL, file="view_me_in_genome_browser", chrom.length.file=NULL) {

	# Variables
	file <- paste(file,".bed", sep="")
	combstates <- model$stateorder
	if (is.null(numstates)) {
		numstates <- length(combstates)
	} else if (numstates > length(combstates)) {
		numstates <- length(combstates)
	}
	if (exclude.state.zero) {
		combstates2use <- setdiff(combstates,0)
		combstates2use <- combstates2use[1:numstates]
		combstates2use <- combstates2use[!is.na(combstates2use)]
		numstates <- length(combstates2use)
	}

	# Collapse the calls
	calls <- cbind(model$coordinates, name=model$states)
	collapsed.calls <- collapse.bins(calls, column2collapseBy=4)
	# Select only desired states
	mask <- rep(FALSE,nrow(collapsed.calls))
	for (istate in combstates2use) {
		mask <- mask | istate==collapsed.calls$name
	}
	collapsed.calls <- collapsed.calls[mask,]

	# Check length of chromosomes if chrom.length.file was given
	if (!is.null(chrom.length.file)) {
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		for (chrom in levels(as.factor(collapsed.calls$chrom))) {
			index <- which(collapsed.calls$chrom == chrom & collapsed.calls$end > chrom.lengths[chrom])
			collapsed.calls$end[index] <- chrom.lengths[chrom]
			if (length(index) >= 1) {
				warning("Adjusted entry in chromosome ",chrom,", because it was longer than the length of the chromosome")
			}
		}
	}

	# Generate the colors for each combinatorial state
	colors <- colors()[grep(colors(), pattern="white|grey|gray|snow", invert=T)]
	step <- length(colors) %/% numstates
	colors <- colors[seq(1,by=step,length=numstates)]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")
	itemRgb <- RGBs[as.factor(collapsed.calls$name)]

	# Write to file
	cat("browser hide all\n", file=file)
	numsegments <- nrow(collapsed.calls)
	df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)

	if (separate.tracks) {
		for (istate in combstates2use) {
			priority <- 100 + which(istate==combstates2use)
			cat(paste("track name=\"comb.state ",istate,"\" description=\"multivariate calls for combinatorial state ",istate,"\" visibility=1 itemRgb=On priority=",priority,"\n", sep=""), file=file, append=TRUE)
			write.table(format(df[df$name==istate,], scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	} else {
		cat(paste("track name=\"comb.state\" description=\"multivariate combinatorial states\" visibility=1 itemRgb=On priority=100\n", sep=""), file=file, append=TRUE)
		write.table(format(df, scientific=FALSE), file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

}

