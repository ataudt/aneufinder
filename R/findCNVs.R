findCNVs <- function(binned.data, ID, eps=0.001, init="random", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, read.cutoff.quantile=0.999) {

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
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")

	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	## Assign variables
# 	state.labels # assigned in global.R
	numstates <- length(state.labels)
	numbins <- length(binned.data)
	reads <- mcols(binned.data)$reads
	iniproc <- which(init==c("standard","random","empiric")) # transform to int

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(reads!=0)) {
		stop("All reads in data are zero. No HMM done.")
	}

	# Filter high reads out, makes HMM faster
	read.cutoff <- as.integer(quantile(reads, read.cutoff.quantile))
	mask <- reads > read.cutoff
	reads[mask] <- read.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		cat(paste0("Replaced read counts > ",read.cutoff," (",names(read.cutoff)," quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'read.cutoff.quantile=1' to disable this filtering.\n"))
	}
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")

		## Initial parameters
		if (init == 'random') {
			A.initial <- matrix(runif(numstates^2), ncol=numstates)
			A.initial <- sweep(A.initial, 1, rowSums(A.initial), "/")			
			proba.initial <- runif(numstates)
			size.initial <- runif(1, min=0, max=1000) * cumsum(ifelse(state.distributions=='dnbinom' | state.distributions=='dbinom', T, F))
			prob.initial <- runif(1) * ifelse(state.distributions=='dnbinom', T, F)
			prob.initial[state.distributions=='dgeom'] <- runif(1)
			lambda.initial <- runif(1, min=0, max=1000) * cumsum(ifelse(state.distributions=='dpois', T, F))
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
		}
	}

	if (num.trials > 1) {

		# Select fit with best loglikelihood
		indexmax <- which.max(unlist(lapply(modellist,"[[","loglik")))
		hmm <- modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		cat("\n\nRerunning try ",indexmax," with eps =",eps,"--------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
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

	# Add useful entries
	hmm$read.cutoff <- read.cutoff
	hmm$ID <- ID
	names(hmm$weights) <- state.labels
	hmm$coordinates <- data.frame(as.character(seqnames(binned.data)), start(ranges(binned.data)), end(ranges(binned.data)))
	names(hmm$coordinates) <- coordinate.names
	hmm$seqlengths <- seqlengths(binned.data)
	class(hmm) <- class.aneufinder.hmm
	hmm$states <- factor(state.labels, levels=state.labels)[hmm$states+1]
	hmm$eps <- eps
	hmm$A <- matrix(hmm$A, ncol=hmm$num.states)
	rownames(hmm$A) <- state.labels
	colnames(hmm$A) <- state.labels
	hmm$distributions <- data.frame()
	hmm$distributions.initial <- data.frame()
	for (idistr in 1:length(hmm$distr.type)) {
		distr <- levels(state.distributions)[hmm$distr.type[idistr]]
		if (distr == 'dnbinom') {
			hmm$distributions <- rbind(hmm$distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], lambda=NA, mu=dnbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dnbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
			hmm$distributions.initial <- rbind(hmm$distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], lambda=NA, mu=dnbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dnbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
		} else if (distr == 'dgeom') {
			hmm$distributions <- rbind(hmm$distributions, data.frame(type=distr, size=NA, prob=hmm$prob[idistr], lambda=NA, mu=dgeom.mean(hmm$prob[idistr]), variance=dgeom.variance(hmm$prob[idistr])))
			hmm$distributions.initial <- rbind(hmm$distributions.initial, data.frame(type=distr, size=NA, prob=hmm$prob.initial[idistr], lambda=NA, mu=dgeom.mean(hmm$prob.initial[idistr]), variance=dgeom.variance(hmm$prob.initial[idistr])))
		} else if (distr == 'delta') {
			hmm$distributions <- rbind(hmm$distributions, data.frame(type=distr, size=NA, prob=NA, lambda=NA, mu=0, variance=0))
			hmm$distributions.initial <- rbind(hmm$distributions.initial, data.frame(type=distr, size=NA, prob=NA, lambda=NA, mu=0, variance=0))
		} else if (distr == 'dpois') {
			hmm$distributions <- rbind(hmm$distributions, data.frame(type=distr, size=NA, prob=NA, lambda=hmm$lambda[idistr], mu=hmm$lambda[idistr], variance=hmm$lambda[idistr]))
			hmm$distributions.initial <- rbind(hmm$distributions.initial, data.frame(type=distr, size=NA, prob=NA, lambda=hmm$lambda.initial[idistr], mu=hmm$lambda.initial[idistr], variance=hmm$lambda.initial[idistr]))
		} else if (distr == 'dbinom') {
			hmm$distributions <- rbind(hmm$distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], lambda=NA, mu=dbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
			hmm$distributions.initial <- rbind(hmm$distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], lambda=NA, mu=dbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
		}
	}
	rownames(hmm$distributions) <- state.labels
	rownames(hmm$distributions.initial) <- state.labels
	hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states)
	rownames(hmm$A.initial) <- state.labels
	colnames(hmm$A.initial) <- state.labels

	# Delete redundant entries
	hmm$size <- NULL
	hmm$prob <- NULL
	hmm$lambda <- NULL
	hmm$size.initial <- NULL
	hmm$prob.initial <- NULL
	hmm$lambda.initial <- NULL
	hmm$use.initial.params <- NULL
	hmm$distr.type <- NULL

	# Issue warnings
	if (num.trials == 1) {
		if (hmm$loglik.delta > hmm$eps) {
			war <- warning("HMM did not converge!\n")
		}
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem. Try again with a lower value for 'read.cutoff.quantile'.")
	} else if (hmm$error == 2) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
	}

	# Return results
	if (!is.null(war)) {
		if (output.if.not.converged == TRUE) {
			return(hmm)
		}
	} else {
		return(hmm)
	}
}
