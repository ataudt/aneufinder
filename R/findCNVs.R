findCNVs <- function(binned.data, ID, method='univariate', eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=10, eps.try=10*eps, num.threads=1, read.cutoff.quantile=0.999, use.gc.corrected.reads=TRUE, strand='*') {

	cat(paste(rep('=',getOption('width')), collapse=''),"\n")
	cat(paste0("Find CNVs for ID = ",ID,"\n"))
	if (method=='univariate') {
		model <- univariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, use.gc.corrected.reads=use.gc.corrected.reads, strand=strand)
	} else if (method=='bivariate') {
		model <- bivariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, use.gc.corrected.reads=use.gc.corrected.reads)
	}

	return(model)

}


univariate.findCNVs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, use.gc.corrected.reads=TRUE, strand='*') {

	## Intercept user input
	IDcheck <- ID  #trigger error if not defined
	if (class(binned.data) != 'GRanges') {
		binned.data <- get(load(binned.data))
		if (class(binned.data) != 'GRanges') stop("argument 'binned.data' expects a GRanges with meta-column 'reads' or a file that contains such an object")
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.logical(use.gc.corrected.reads)!=0) stop("argument 'use.gc.corrected.reads' expects a logical (TRUE or FALSE)")
	if (check.strand(strand)!=0) stop("argument 'strand' expects either '+', '-' or '*'")

	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	## Assign variables
# 	state.labels # assigned in global.R
	numstates <- length(state.labels)
	numbins <- length(binned.data)
	iniproc <- which(init==c("standard","random")) # transform to int
	if (strand=='+') {
		select <- 'preads'
	} else if (strand=='-') {
		select <- 'mreads'
	} else if (strand=='*') {
		select <- 'reads'
	}
	if (use.gc.corrected.reads) {
		if (paste0(select,'.gc') %in% names(mcols(binned.data))) {
			select <- paste0(select,'.gc')
		} else {
			warning("Cannot use GC-corrected reads because they are not in the binned data. Continuing without GC-correction.")
		}
	}
	reads <- mcols(binned.data)[,select]

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(reads!=0)) {
		stop("All reads in data are zero. No HMM done.")
	}

	# Filter high reads out, makes HMM faster
	read.cutoff <- quantile(reads, read.cutoff.quantile)
	names.read.cutoff <- names(read.cutoff)
	read.cutoff <- as.integer(read.cutoff)
	mask <- reads > read.cutoff
	reads[mask] <- read.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		cat(paste0("Replaced read counts > ",read.cutoff," (",names.read.cutoff," quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'read.cutoff.quantile=1' to disable this filtering. This filtering was done to increase the speed of the HMM and should not affect the results.\n\n"))
	}
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat(paste0("Trial ",i_try," / ",num.trials,"\n"))

		## Initial parameters
		if (init == 'random') {
			A.initial <- matrix(runif(numstates^2), ncol=numstates)
			A.initial <- sweep(A.initial, 1, rowSums(A.initial), "/")			
			proba.initial <- runif(numstates)
			size.initial <- runif(1, min=0, max=100) * cumsum(ifelse(state.distributions=='dnbinom' | state.distributions=='dbinom', T, F))
			prob.initial <- runif(1) * ifelse(state.distributions=='dnbinom', T, F)
			prob.initial[state.distributions=='dgeom'] <- runif(1)
			lambda.initial <- runif(1, min=0, max=100) * cumsum(ifelse(state.distributions=='dpois', T, F))
		} else if (init == 'standard') {
			A.initial <- matrix(NA, ncol=numstates, nrow=numstates)
			for (irow in 1:numstates) {
				for (icol in 1:numstates) {
					if (irow==icol) { A.initial[irow,icol] <- 0.9 }
					else { A.initial[irow,icol] <- 0.1/(numstates-1) }
				}
			}
			proba.initial <- rep(1/numstates, numstates)
			mean.initial <- mean(reads[reads>0])/2 * cumsum(ifelse(state.distributions=='dnbinom' | state.distributions=='dbinom', T, F))
			var.initial <- var(reads[reads>0])/2 * cumsum(ifelse(state.distributions=='dnbinom' | state.distributions=='dbinom', T, F))
			size.initial <- rep(0,numstates)
			prob.initial <- rep(0,numstates)
			mask <- state.distributions=='dnbinom'
			size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
			prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
			mask <- state.distributions=='dbinom'
			size.initial[mask] <- dbinom.size(mean.initial[mask], var.initial[mask])
			prob.initial[mask] <- dbinom.prob(mean.initial[mask], var.initial[mask])
			prob.initial[state.distributions=='dgeom'] <- 0.9
			lambda.initial <- mean(reads[reads>0])/2 * cumsum(ifelse(state.distributions=='dpois', T, F))
		}
	
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			state.labels = as.integer(state.labels), # int* state_labels
			size = double(length=numstates), # double* size
			prob = double(length=numstates), # double* prob
			lambda = double(length=numstates), # double* lambda
			num.iterations = as.integer(max.iter), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps.try), # double* eps
			states = integer(length=numbins), # int* states
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			distr.type = as.integer(state.distributions), # int* distr_type
			ini.proc = as.integer(iniproc), # int* iniproc
			size.initial = as.vector(size.initial), # double* initial_size
			prob.initial = as.vector(prob.initial), # double* initial_prob
			lambda.initial = as.vector(lambda.initial), # double* initial_lambda
			A.initial = as.vector(A.initial), # double* initial_A
			proba.initial = as.vector(proba.initial), # double* initial_proba
			use.initial.params = as.logical(1), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)

		hmm$eps <- eps.try
		if (num.trials > 1) {
			if (hmm$loglik.delta > hmm$eps) {
				warning("HMM did not converge in trial run ",i_try,"!\n")
			}
			# Store model in list
			hmm$reads <- NULL
			modellist[[i_try]] <- hmm
			# Check if disomic state is most frequent
			idisomy <- which(state.labels=='disomy')
			if (which.max(hmm$weights)==idisomy) {
				break
			}
			init <- 'random'
		}
	}

	if (num.trials > 1) {

		# Select fit with highest weight in state disomic
		# Mathematically we should select the fit with highest loglikelihood. If we think the fit with the highest loglikelihood is incorrect, we should change the underlying model. However, this is very complex and we choose to select a fit that we think is (more) correct, although it has not the highest support given our (imperfect) model.
		idisomy <- which(state.labels=='disomy')
		indexmax <- which.max(unlist(lapply(lapply(modellist,'[[','weights'), '[', idisomy)))
		hmm <- modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		cat(paste0("Rerunning trial ",indexmax," with eps = ",eps,"\n"))
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			state.labels = as.integer(state.labels), # int* state_labels
			size = double(length=numstates), # double* size
			prob = double(length=numstates), # double* prob
			lambda = double(length=numstates), # double* lambda
			num.iterations = as.integer(max.iter), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps), # double* eps
			states = integer(length=numbins), # int* states
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			distr.type = as.integer(state.distributions), # int* distr_type
			ini.proc = as.integer(iniproc), # int* iniproc
			size.initial = as.vector(hmm$size), # double* initial_size
			prob.initial = as.vector(hmm$prob), # double* initial_prob
			lambda.initial = as.vector(hmm$lambda), # double* initial_lambda
			A.initial = as.vector(hmm$A), # double* initial_A
			proba.initial = as.vector(hmm$proba), # double* initial_proba
			use.initial.params = as.logical(1), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)
	}

	### Make return object ###
		result <- list()
		result$ID <- ID
	## Bin coordinates and states ###
		result$bins <- binned.data
		result$bins$state <- state.labels[hmm$states]
	## Segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- result$bins
		red.gr.list <- GRangesList()
		for (state in state.labels) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		result$segments <- red.gr
		seqlengths(result$segments) <- seqlengths(binned.data)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
	## Parameters
		# Weights
		result$weights <- hmm$weights
		names(result$weights) <- state.labels
		# Transition matrices
		transitionProbs <- matrix(hmm$A, ncol=hmm$num.states)
		rownames(transitionProbs) <- state.labels
		colnames(transitionProbs) <- state.labels
		result$transitionProbs <- transitionProbs
		transitionProbs.initial <- matrix(hmm$A.initial, ncol=hmm$num.states)
		rownames(transitionProbs.initial) <- state.labels
		colnames(transitionProbs.initial) <- state.labels
		result$transitionProbs.initial <- transitionProbs.initial
		# Initial probs
		result$startProbs <- hmm$proba
		names(result$startProbs) <- paste0("P(",state.labels,")")
		result$startProbs.initial <- hmm$proba.initial
		names(result$startProbs.initial) <- paste0("P(",state.labels,")")
		# Distributions
			distributions <- data.frame()
			distributions.initial <- data.frame()
			for (idistr in 1:length(hmm$distr.type)) {
				distr <- levels(state.distributions)[hmm$distr.type[idistr]]
				if (distr == 'dnbinom') {
					distributions <- rbind(distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], lambda=NA, mu=dnbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dnbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
					distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], lambda=NA, mu=dnbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dnbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
				} else if (distr == 'dgeom') {
					distributions <- rbind(distributions, data.frame(type=distr, size=NA, prob=hmm$prob[idistr], lambda=NA, mu=dgeom.mean(hmm$prob[idistr]), variance=dgeom.variance(hmm$prob[idistr])))
					distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=NA, prob=hmm$prob.initial[idistr], lambda=NA, mu=dgeom.mean(hmm$prob.initial[idistr]), variance=dgeom.variance(hmm$prob.initial[idistr])))
				} else if (distr == 'delta') {
					distributions <- rbind(distributions, data.frame(type=distr, size=NA, prob=NA, lambda=NA, mu=0, variance=0))
					distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=NA, prob=NA, lambda=NA, mu=0, variance=0))
				} else if (distr == 'dpois') {
					distributions <- rbind(distributions, data.frame(type=distr, size=NA, prob=NA, lambda=hmm$lambda[idistr], mu=hmm$lambda[idistr], variance=hmm$lambda[idistr]))
					distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=NA, prob=NA, lambda=hmm$lambda.initial[idistr], mu=hmm$lambda.initial[idistr], variance=hmm$lambda.initial[idistr]))
				} else if (distr == 'dbinom') {
					distributions <- rbind(distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], lambda=NA, mu=dbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
					distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], lambda=NA, mu=dbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
				}
			}
			rownames(distributions) <- state.labels
			rownames(distributions.initial) <- state.labels
			result$distributions <- distributions
			result$distributions.initial <- distributions
	## Convergence info
		convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec)
		result$convergenceInfo <- convergenceInfo
	## Add class
		class(result) <- class.univariate.hmm

	### Issue warnings ###
	if (num.trials == 1) {
		if (hmm$loglik.delta > eps) {
			war <- warning("HMM did not converge!\n")
		}
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem. Try again with a lower value for 'read.cutoff.quantile'.")
	} else if (hmm$error == 2) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
	}

	# Return results
	return(result)
}


bivariate.findCNVs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, use.gc.corrected.reads=TRUE) {

	## Intercept user input
	IDcheck <- ID  #trigger error if not defined
	if (class(binned.data) != 'GRanges') {
		binned.data <- get(load(binned.data))
		if (class(binned.data) != 'GRanges') stop("argument 'binned.data' expects a GRanges with meta-column 'reads' or a file that contains such an object")
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.logical(use.gc.corrected.reads)!=0) stop("argument 'use.gc.corrected.reads' expects a logical (TRUE or FALSE)")

	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	### Split into strands and run univariate findCNVs
	cat('\n')
	cat(paste(rep('-',getOption('width')), collapse=''),"\n")
	cat(paste0("Running '+'-strand\n"))
	plus.model <- univariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, use.gc.corrected.reads=use.gc.corrected.reads, strand='+')
	cat('\n')
	cat(paste(rep('-',getOption('width')), collapse=''),"\n")
	cat(paste0("Running '-'-strand\n"))
	minus.model <- univariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, use.gc.corrected.reads=use.gc.corrected.reads, strand='-')
	models <- list(plus=plus.model, minus=minus.model)

	### Prepare the multivariate HMM
	cat('\n')
	cat(paste(rep('-',getOption('width')), collapse=''),"\n")
	cat(paste0("Preparing bivariate HMM\n\n"))

		## Extract reads and other stuff
		distributions <- lapply(models, '[[', 'distributions')
		weights <- lapply(models, '[[', 'weights')
		num.uni.states <- length(models[[1]]$weights)
		num.models <- length(models)
		num.bins <- length(models[[1]]$bins)
		comb.states <- paste(levels(models[[1]]$bins$state), levels(models[[2]]$bins$state))
		comb.states <- factor(comb.states, levels=comb.states)
		states.list <- list()
		for (model in models) { states.list[[length(states.list)+1]] <- model$bins$state }
		comb.states.per.bin <- factor(do.call(paste, states.list), levels=levels(comb.states))
		num.comb.states <- length(comb.states)
		select <- 'reads'
		if (use.gc.corrected.reads) {
			if (paste0(select,'.gc') %in% names(mcols(binned.data))) {
				select <- paste0(select,'.gc')
			} else {
				warning("Cannot use GC-corrected reads because they are not in the binned data. Continuing without GC-correction.")
			}
		}
		reads <- matrix(c(mcols(models$plus$bins)[,paste0('p',select)], mcols(models$minus$bins)[,paste0('m',select)]), ncol=num.models, dimnames=list(bin=1:num.bins, strand=names(models)))
		maxreads <- max(reads)

		## Pre-compute z-values for each number of reads
		cat("Computing pre z-matrix...")
		z.per.read <- array(NA, dim=c(maxreads+1, num.models, num.uni.states))
		xreads <- 0:maxreads
		for (istrand in 1:num.models) {
			model <- models[[istrand]]
			for (istate in 1:num.uni.states) {
				if (model$distributions[istate,'type']=='dnbinom') {
					size <- model$distributions[istate,'size']
					prob <- model$distributions[istate,'prob']
					u <- pnbinom(xreads, size, prob)
				} else if (model$distributions[istate,'type']=='delta') {
					u <- rep(1, length(xreads))
				}
				qnorm_u <- qnorm(u)
				mask <- qnorm_u==Inf
				qnorm_u[mask] <- qnorm(1-1e-16)
				z.per.read[ , istrand, istate] <- qnorm_u
			}
		}
		cat(" done\n")

		## Compute the z matrix
		cat("Transfering values into z-matrix...")
		z.per.bin = array(NA, dim=c(num.bins, num.models, num.uni.states), dimnames=list(bin=1:num.bins, strand=names(models), levels(models[[1]]$bins$state)))
		for (istrand in 1:num.models) {
			model <- models[[istrand]]
			for (istate in 1:num.uni.states) {
				z.per.bin[ , istrand, istate] = z.per.read[reads[,istrand]+1, istrand, istate]
			}
		}
		remove(z.per.read)
		cat(" done\n")

		## Calculate correlation matrix
		cat("Computing inverse of correlation matrix...")
		correlationMatrix <- array(0, dim=c(num.models,num.models,num.comb.states), dimnames=list(strand=names(models), strand=names(models), comb.state=comb.states))
		correlationMatrixInverse <- correlationMatrix
		determinant <- rep(0, num.comb.states)
		names(determinant) <- comb.states
		usestateTF <- rep(NA, num.comb.states) # TRUE, FALSE vector for usable states
		names(usestateTF) <- comb.states
		for (state in comb.states) {
			istate <- strsplit(state, ' ')[[1]]
			mask <- which(comb.states.per.bin==state)
			# Subselect z
			z.temp <- matrix(NA, ncol=length(istate), nrow=length(mask))
			for (i1 in 1:length(istate)) {
				z.temp[,i1] <- z.per.bin[mask, i1, istate[i1]]
			}
			temp <- tryCatch({
				correlationMatrix[,,state] <- cor(z.temp)
				determinant[state] <- det( correlationMatrix[,,state] )
				correlationMatrixInverse[,,state] <- solve(correlationMatrix[,,state])
				if (length(z.temp)>0) {
					usestateTF[state] <- TRUE
				} else {
					usestateTF[state] <- FALSE
				}
			}, warning = function(war) {
				usestateTF[state] <<- FALSE
				war
			}, error = function(err) {
				usestateTF[state] <<- FALSE
				err
			})
		}
		cat(" done\n")

		# Use nullsomy state anyways
		usestateTF[1] <- TRUE
		correlationMatrix = correlationMatrix[,,usestateTF]
		correlationMatrixInverse = correlationMatrixInverse[,,usestateTF]
		comb.states = comb.states[usestateTF]
		comb.states <- droplevels(comb.states)
		determinant = determinant[usestateTF]
		num.comb.states <- length(comb.states)

		## Calculate multivariate densities for each state
		cat("Calculating multivariate densities...")
		densities <- matrix(1, ncol=num.comb.states, nrow=num.bins, dimnames=list(bin=1:num.bins, comb.state=comb.states))
		for (state in comb.states) {
			i <- which(state==comb.states)
			istate <- strsplit(state, ' ')[[1]]
			product <- 1
			for (istrand in 1:num.models) {
				if (models[[istrand]]$distributions[istate[istrand],'type'] == 'dnbinom') {
					exponent <- -0.5 * apply( ( z.per.bin[ , , i] %*% (correlationMatrixInverse[ , , i] - diag(num.models)) ) * z.per.bin[ , , i], 1, sum)
					size <- models[[istrand]]$distributions[istate[istrand],'size']
					prob <- models[[istrand]]$distributions[istate[istrand],'prob']
					product <- product * dnbinom(reads[,istrand], size, prob)
					densities[,i] <- product * determinant[i]^(-0.5) * exp( exponent )
				}
				else if (models[[istrand]]$distributions[istate[istrand],'type'] == 'delta') {
					densities[,i] <- densities[,i] * ifelse(reads[,istrand]==0, 1, 0)
				}
			}
		}
		# Check if densities are > 1
# 		if (any(densities>1)) stop("Densities > 1")
# 		if (any(densities<0)) stop("Densities < 0")
		densities[densities>1] <- 1
		densities[densities<0] <- 0
		# Check if densities are 0 everywhere in some bins
		check <- which(apply(densities, 1, sum) == 0)
		if (length(check)>0) {
			if (check[1]==1) {
				densities[1,] <- rep(1e-10, ncol(densities))
				check <- check[-1]
			}
			for (icheck in check) {
				densities[icheck,] <- densities[icheck-1,]
			}
		}
		cat(" done\n")
		

	### Run the multivariate HMM
	# Call the C function
	use.initial <- FALSE
	hmm <- .C("R_multivariate_hmm",
		densities = as.double(densities), # double* D
		num.bins = as.integer(num.bins), # int* T
		num.comb.states = as.integer(num.comb.states), # int* N
		num.strands = as.integer(num.models), # int* Nmod
		comb.states = as.integer(comb.states), # int* comb_states
		num.iterations = as.integer(max.iter), # int* maxiter
		time.sec = as.integer(max.time), # double* maxtime
		loglik.delta = as.double(eps), # double* eps
		states = integer(length=num.bins), # int* states
		A = double(length=num.comb.states*num.comb.states), # double* A
		proba = double(length=num.comb.states), # double* proba
		loglik = double(length=1), # double* loglik
		A.initial = double(length=num.comb.states*num.comb.states), # double* initial_A
		proba.initial = double(length=num.comb.states), # double* initial_proba
		use.initial.params = as.logical(use.initial), # bool* use_initial_params
		num.threads = as.integer(num.threads), # int* num_threads
		error = as.integer(0) # error handling
		)
			
	### Check convergence ###
	war <- NULL
	if (hmm$loglik.delta > eps) {
		war <- warning("HMM did not converge!\n")
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
	} else if (hmm$error == 2) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
	}

	### Make return object ###
		result <- list()
		result$ID <- ID
	## Bin coordinates and states
		result$bins <- binned.data
		result$bins$state <- comb.states[hmm$states]
		# Get states as factors in data.frame
		matrix.states <- matrix(unlist(strsplit(as.character(result$bins$state), split=' ')), byrow=T, ncol=num.models, dimnames=list(bin=1:num.bins, strand=names(models)))
		names <- c('pstate','mstate')
		for (i1 in 1:num.models) {
			mcols(result$bins)[names[i1]] <- factor(matrix.states[,i1], levels=levels(models[[i1]]$bins$state))
		}
	## Segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- result$bins
		red.gr.list <- GRangesList()
		for (state in comb.states) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		result$segments <- red.gr
		seqlengths(result$segments) <- seqlengths(result$bins)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		## Parameters
			# Weights
			tstates <- table(result$bins$state)
			result$weights <- tstates/sum(tstates)
			# Transition matrices
			result$transitionProbs <- matrix(hmm$A, ncol=num.comb.states)
			colnames(result$transitionProbs) <- comb.states
			rownames(result$transitionProbs) <- comb.states
			result$transitionProbs.initial <- matrix(hmm$A.initial, ncol=num.comb.states)
			colnames(result$transitionProbs.initial) <- comb.states
			rownames(result$transitionProbs.initial) <- comb.states
			# Initial probs
			result$startProbs <- hmm$proba
			names(result$startProbs) <- paste0("P(",comb.states,")")
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- paste0("P(",comb.states,")")
			# Distributions
			result$distributions <- distributions
			names(result$distributions) <- names(models)
		## Convergence info
			convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec)
			result$convergenceInfo <- convergenceInfo
		## Correlation matrices
# 			result$correlation.matrix <- correlationMatrix2use
		## Add class
			class(result) <- class.multivariate.hmm

	return(result)

}


