# ==================================================================================
# Generate a data.frame with majority-state for each chromosome and input file/model
# ==================================================================================
get.state.table <- function(hmm.list, numCPU=1) {

	## Transform to GRanges
# 	grlist <- lapply(hmm.list, hmm2GRanges, reduce=F)
	grlist <- hmmList2GRangesList(hmm.list, reduce=F, numCPU=numCPU)

	## Split by chromosome
	grsplitlist <- lapply(grlist, function(gr) { split(gr, seqnames(gr)) })

	## Get the tables per model/file[[]] and chromosome[[]][[]]
	tablell <- lapply(grsplitlist, function(grsplit) { lapply(grsplit, function(gr) { table(mcols(gr)$state) }) })

	## Transform to array. Dimensions are [model/file, chromosome, state]
	states <- array(dim=c(length(hmm.list), length(tablell[[1]]), length(tablell[[1]][[1]])), dimnames=list('model'=1:length(tablell), 'chromosome'=names(tablell[[1]]), 'state'=names(tablell[[1]][[1]])))
	for (i1 in 1:dim(states)[1]) {
		for (i2 in 1:dim(states)[2]) {
			tc <- tryCatch({
				states[i1,i2,] <- tablell[[i1]][[i2]]
			}, error = function(err) {
				err
			})
		}
	}

	## Get the majority state per chromosome of each model/file
	majorstate.f <- apply(states, c(1,2), function(x) { names(x)[which.max(x)] })
	majorstate <- matrix(NA, ncol=ncol(majorstate.f), nrow=nrow(majorstate.f))
	colnames(majorstate) <- colnames(majorstate.f)
	rownames(majorstate) <- rownames(majorstate.f)
	for (i1 in 1:ncol(majorstate.f)) {
		majorstate[,i1] <- factor(majorstate.f[,i1], levels=state.labels)
	}
	majorstate.f <- matrix(state.labels[majorstate], ncol=ncol(majorstate), nrow=nrow(majorstate))
	colnames(majorstate.f) <- colnames(majorstate)
	rownames(majorstate.f) <- rownames(majorstate)

	## Plot the percentages of aneuploidies/states per chromosome
	df <- melt(majorstate.f, varnames=c("model","chromosome"), value.name="state")
	df$state <- factor(df$state, levels=state.labels)
	ggplt <- ggplot(df) + geom_bar(aes(x=chromosome, fill=state, y=..count..)) + theme_bw() + ylab('number of samples') + scale_fill_manual(values=state.colors)

	## Return results
	l <- list(table=majorstate, plot=ggplt)
	return(l)


}


# ==================================================================================
# Build a dendrogram
# ==================================================================================
get.dendrogram <- function(hmm.list, numCPU=1) {

	### Generate the GRanges consensus template with variable bins that can hold the states of each model
	## Load the files
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges in reduced representation
	temp <- hmmList2GRangesList(hmm.list, reduce=TRUE, numCPU=numCPU, consensus=TRUE)
	grlred <- temp$grl
	consensus <- temp$consensus

	## Split into non-overlapping fragments
	## Overlap each models' states with that consensus template
	cat('calculate overlap\n')
	constates <- foreach (gr = grlred, .packages='GenomicRanges', .combine='cbind') %do% {
		splt <- split(gr, mcols(gr)$state)
		mind <- as.matrix(findOverlaps(consensus, splt))
		col <- matrix(-1, nrow=length(consensus), ncol=1)
		col[mind[,'queryHits'],1] <- mind[,'subjectHits']
		col
	}
	colnames(constates) <- unlist(lapply(hmm.list, '[[', 'ID'))
		
	## Distance measure
	cat('calculating distance\n')
	wcor <- cov.wt(constates, wt=as.numeric(width(consensus)), cor=T)
	dist <- as.dist(1-wcor$cor)
	## Dendrogram
	hc <- hclust(dist)

	## Return results
	return(hc)

}


# ==================================================================================
# Get a simple measure of reproducibility between a multiple (num.query.samples)
# cell sample (query) and single cell (subjects) samples
# ==================================================================================
get.reproducibility.one2many <- function(query, num.query.samples=1, subjects, num.tests=10, numCPU=1) {

	if (check.univariate.model(query)!=0) {
		cat("loading univariate HMM from file\n")
		query <- get(load(query))
		if (check.univariate.model(query)!=0) stop("argument 'query' expects a univariate hmm or a file that contains a univariate hmm")
	}

	## Transform to GRanges
	query.gr <- hmm2GRanges(query, reduce=F)
	subjects.grl <- hmmList2GRangesList(subjects, reduce=T, numCPU=numCPU)

	## Overlap each subject's states with query
	cat('overlapping each subject\'s states with query\n')
	cl <- makeCluster(numCPU)
	registerDoParallel(cl)
	constates <- foreach (subject.gr = subjects.grl, .packages='GenomicRanges', .combine='cbind') %dopar% {
		splt <- split(subject.gr, mcols(subject.gr)$state)
		mind <- as.matrix(findOverlaps(query.gr, splt, select='first'))
	}
	stopCluster(cl)

	## Comparing sample of subjects to query
	cat('comparing samples of subjects to query\n')
	means <- matrix(NA, nrow=num.tests, ncol=length(levels(mcols(query.gr)$state)))
	colnames(means) <- levels(mcols(query.gr)$state)
	meanconstates <- matrix(NA, nrow=nrow(constates), ncol=num.tests)
	for (itest in 1:num.tests) {
		cat('test',itest,'\n')

		## Sample the columns to compare to the query
		cols <- sample(1:length(subjects), num.query.samples)

		## Mean of each position
		meanconstates[,itest] <- apply(constates[,cols], 1, mean, na.rm=T)

		## Subjects-mean of each query-level
		meanlevel <- NULL
		varlevel <- NULL
		for (level in levels(mcols(query.gr)$state)) {
			meanlevel[level] <- mean(meanconstates[mcols(query.gr)$state==level])
			varlevel[level] <- var(meanconstates[mcols(query.gr)$state==level])
		}
		means[itest,] <- meanlevel
	}
	## Plot the results
	df <- melt(data.frame(state=mcols(query.gr)$state, meanstate=meanconstates))
	ggplt <- ggplot(df) + geom_boxplot(aes(x=state, y=value, fill=variable)) + theme_bw() + coord_cartesian(ylim=c(1,length(levels(mcols(query.gr)$state)))) + scale_y_continuous(labels=levels(mcols(query.gr)$state)) + xlab(paste0(num.query.samples,' cell sample')) + ylab(paste0('single cell samples'))

	df <- melt(means, varnames=c('test','state'))
	ggplt1 <- ggplot(df) + geom_boxplot(aes(x=state, y=value)) + theme_bw() + coord_cartesian(ylim=c(1,length(levels(mcols(query.gr)$state)))) + scale_y_continuous(labels=levels(mcols(query.gr)$state)) + xlab(paste0(num.query.samples,' cell sample')) + ylab(paste0('single cell samples'))

	## Comparing all subjects to query (no means)
	cat('comparing all samples to query\n')
	df.constates <- as.data.frame(cbind(mcols(query.gr)$state, constates))
	names(df.constates)[1] <- 'query'
	df <- melt(df.constates, id.vars='query', value.name='state')
	df <- df[seq(from=1, to=nrow(df), by=10000),]
	df$query <- factor(levels(mcols(query.gr)$state)[df$query], levels=levels(mcols(query.gr)$state))

	ggplt2 <- ggplot(df) + geom_boxplot(aes(x=query, y=state)) + theme_bw() + coord_cartesian(ylim=c(1,length(levels(mcols(query.gr)$state)))) + scale_y_continuous(labels=levels(mcols(query.gr)$state)) + xlab(paste0(num.query.samples,' cell sample')) + ylab(paste0('single cell samples'))

	df$state <- factor(levels(mcols(query.gr)$state)[df$state], levels=levels(mcols(query.gr)$state))
	ggplt3 <- ggplot(df) + geom_jitter(aes(x=query, y=state), position=position_jitter(width=0.3, height=0.2)) + theme_bw() + coord_cartesian(ylim=c(1,length(levels(mcols(query.gr)$state)))) + xlab(paste0(num.query.samples,' cell sample')) + ylab(paste0('single cell samples'))
	
	
	


	## Return results
	l <- list(meanplot=ggplt, meanmeanplot=ggplt1, allbox=ggplt2, alljitter=ggplt3)
	return(l)



}



# ==================================================================================
# Get a simple measure of reproducibility between sets of single cell samples
# ==================================================================================
get.reproducibility.many2many <- function(samples, num.tests=10, size.test=10, numCPU=1) {

	## Transform to GRanges
	temp <- hmmList2GRangesList(samples, numCPU=numCPU, reduce=T, consensus=T)
	samples.grl <- temp$grl
	consensus <- temp$consensus
	constates <- temp$constates

	## Comparing samples against each other
	cat('comparing samples against each other\n')
	meanconstates <- array(dim=c(length(consensus), 2, num.tests))
	for (itest in 1:num.tests) {
		cat('test',itest,'\n')

		## Sample the columns to compare to the query
		cols <- sample(1:length(samples), 2*size.test)
		cols1 <- cols[1:size.test]
		cols2 <- cols[(size.test+1):(2*size.test)]

		## Mean of each position
		if (size.test==1) {
			meanconstates[,1,itest] <- constates[,cols1]
			meanconstates[,2,itest] <- constates[,cols2]
		} else {
			meanconstates[,1,itest] <- apply(constates[,cols1], 1, mean, na.rm=T)
			meanconstates[,2,itest] <- apply(constates[,cols2], 1, mean, na.rm=T)
		}

	}
	## Mean along all tests
	means <- apply(meanconstates, c(1,2), mean)
	means <- cbind(as.data.frame(means), width=width(consensus))

	## Plot the results
	states <- levels(mcols(samples.grl[[1]])$state)
	ggplt1 <- ggplot(as.data.frame(means)) + geom_point(aes(x=V1, y=V2, alpha=width)) + theme_bw() + coord_cartesian(ylim=c(1,length(states))) + scale_x_continuous(breaks=1:length(states), labels=states) + scale_y_continuous(breaks=1:length(states), labels=states) + xlab(paste0(size.test,' single cell samples')) + ylab(paste0(size.test,' single cell samples'))

	## Return results
	return(ggplt1)

}
