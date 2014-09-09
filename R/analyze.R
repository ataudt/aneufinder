# ==================================================================================
# Generate a data.frame with majority-state for each chromosome and input file/model
# ==================================================================================
get.state.table <- function(modellist) {

	## Intercept user input
	if (check.univariate.modellist(modellist)!=0) {
		cat("Loading univariate HMMs from files ...")
		mlist <- NULL
		for (modelfile in modellist) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		modellist <- mlist
		remove(mlist)
		cat(" done\n")
		if (check.univariate.modellist(modellist)!=0) stop("argument 'modellist' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}
	
	# -----------------------------------
	# Percentage of aneuploid chromosomes
	# -----------------------------------

	## Transform to GRanges
	grlist <- lapply(modellist, hmm2GRanges, reduce=F)

	## Split by chromosome
	grsplitlist <- lapply(grlist, function(gr) { split(gr, seqnames(gr)) })

	## Get the tables per model/file[[]] and chromosome[[]][[]]
	tablell <- lapply(grsplitlist, function(grsplit) { lapply(grsplit, function(gr) { table(mcols(gr)$state) }) })

	## Transform to array. Dimensions are [model/file, chromosome, state]
	states <- array(dim=c(length(modellist), length(tablell[[1]]), length(tablell[[1]][[1]])), dimnames=list('model'=1:length(tablell), 'chromosome'=names(tablell[[1]]), 'state'=names(tablell[[1]][[1]])))
	for (i1 in 1:dim(states)[1]) {
		for (i2 in 1:dim(states)[2]) {
			states[i1,i2,] <- tablell[[i1]][[i2]]
		}
	}

	## Get the majority state per chromosome of each model/file
	majorstate.f <- apply(states, c(1,2), function(x) { names(x)[which.max(x)] })
	factor2number <- 1:length(state.labels) -1
	names(factor2number) <- state.labels
	majorstate <- majorstate.f
	for (i1 in 1:ncol(majorstate.f)) {
		majorstate[,i1] <- factor2number[as.character(majorstate.f[,i1])]
	}
	majorstate <- as.data.frame(majorstate)

	## Plot the percentages of aneuploidies/states per chromosome
	library(ggplot2)
	library(reshape2)
	df <- melt(majorstate.f, measure.vars=colnames(majorstate), variable.name="chromosome", value.name="state")
	df$state <- factor(df$state, levels=state.labels)
	ggplt <- ggplot(df) + geom_bar(aes(n=nrow(majorstate), x=chromosome, fill=state, y=..count../n)) + theme_bw() + ylab('number of samples') + scale_fill_manual(values=state.colors) + scale_y_continuous(labels=scales::percent_format())

	# ----------
	# Dendrogram
	# ----------

	### Generate the GRanges consensus template with variable bins that can hold the states of each model
	## Transform to GRanges in reduced representation
	grlred <- GRangesList()
	for (i1 in 1:length(modellist)) {
		grlred[[i1]] <- hmm2GRanges(modellist[[i1]], reduce=T)
	}
	## Split into non-overlapping fragments
	consensus <- disjoin(unlist(grlred))
	## Overlap each models' states with that consensus template
	constates <- matrix(NA, ncol=length(modellist), nrow=length(consensus))
	for (i1 in 1:length(modellist)) {
		splt <- split(grlred[[i1]], mcols(grlred[[i1]])$state)
		mind <- as.matrix(findOverlaps(consensus, splt))
		constates[,i1] <- mind[,'subjectHits']
	}
	## Distance measure
	wcor <- cov.wt(constates, wt=as.numeric(width(consensus)), cor=T)
	dist <- as.dist(1-wcor$cor)
	## Dendrogram
	hc <- hclust(dist)

	## Return results
	l <- list(table=majorstate, percentage.plot=ggplt, dendrogram=hc)
	return(l)


}

