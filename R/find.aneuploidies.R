find.aneuploidies <- function(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE) {

	## Intercept user input
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.it)!=0) stop("argument 'max.it' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")


	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	names(binned.data) <- binned.data.names # defined globally outside this function

	## Assign variables
# 	state.labels # assigned globally outside this function
	numstates <- length(state.labels)
	numbins <- length(binned.data$reads)

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(binned.data$reads!=0)) {
		stop("All reads in data are zero. No univariate HMM done.")
	}

	# Filter high reads out, makes HMM faster
	if (filter.reads) {
		limit <- 10*ceiling(var(binned.data$reads))
		mask <- binned.data$reads > limit
		binned.data$reads[mask] <- limit
		numfiltered <- length(which(mask))
		if (numfiltered > 0) {
			warning(paste("There are very high read counts (probably artificial) in your data. Replaced read counts > ",limit," (10*variance) by ",limit," in ",numfiltered," bins. Set option 'filter.reads=FALSE' to disable this filtering.", sep=""))
		}
	}
	
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(binned.data$reads), # double* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			means = double(length=numstates), # double* means
			variances = double(length=numstates), # double* variances
			num.iterations = as.integer(max.it), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps.try), # double* eps
			posteriors = double(length=numbins * numstates), # double* posteriors
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			means.initial = double(length=numstates), # double* initial_means
			variances.initial = double(length=numstates), # double* initial_variances
			A.initial = double(length=numstates*numstates), # double* initial_A
			proba.initial = double(length=numstates), # double* initial_proba
			use.initial.params = as.logical(0), # bool* use_initial_params
			num.threads = as.integer(num.threads) # int* num_threads
		)

		hmm$eps <- eps.try
		hmm$A <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
		rownames(hmm$A) <- state.labels
		colnames(hmm$A) <- state.labels
		hmm$distributions <- cbind(mean=hmm$means, variance=hmm$variances, size=fsize(hmm$means,hmm$variances), prob=fprob(hmm$means,hmm$variances))
		rownames(hmm$distributions) <- state.labels
		hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
		rownames(hmm$A.initial) <- state.labels
		colnames(hmm$A.initial) <- state.labels
		hmm$distributions.initial <- cbind(mean=hmm$means.initial, variance=hmm$variances.initial, size=fsize(hmm$means.initial,hmm$variances.initial), prob=fprob(hmm$means.initial,hmm$variances.initial))
		rownames(hmm$distributions.initial) <- state.labels
		if (num.trials > 1) {
			if (hmm$loglik.delta > hmm$eps) {
				warning("HMM did not converge in trial run ",i_try,"!\n")
			}
			# Store model in list
			hmm$posteriors <- NULL
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
			reads <- as.integer(binned.data$reads), # double* O
			num.bins <- as.integer(numbins), # int* T
			num.states <- as.integer(numstates), # int* N
			means <- double(length=numstates), # double* means
			variances <- double(length=numstates), # double* variances
			num.iterations <- as.integer(max.it), #  int* maxiter
			time.sec <- as.integer(max.time), # double* maxtime
			loglik.delta <- as.double(eps), # double* eps
			posteriors <- double(length=numbins * numstates), # double* posteriors
			A <- double(length=numstates*numstates), # double* A
			proba <- double(length=numstates), # double* proba
			loglik <- double(length=1), # double* loglik
			weights <- double(length=numstates), # double* weights
			means.initial <- as.vector(hmm$distributions[,'mean']), # double* initial_means
			variances.initial <- as.vector(hmm$distributions[,'variance']), # double* initial_variances
			A.initial <- as.vector(hmm$A), # double* initial_A
			proba.initial <- as.vector(hmm$proba), # double* initial_proba
			use.initial.params <- as.logical(1), # bool* use_initial_params
			num.threads <- as.integer(num.threads) # int* num_threads
		)
	}

	# Add useful entries
	hmm$coordinates <- binned.data[,coordinate.names]
	hmm$posteriors <- matrix(hmm$posteriors, ncol=hmm$num.states)
	colnames(hmm$posteriors) <- paste("P(",state.labels,")", sep="")
	class(hmm) <- "aneufinder.model"
	hmm$states <- state.labels[apply(hmm$posteriors, 1, which.max)]
	hmm$eps <- eps
	hmm$A <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A) <- state.labels
	colnames(hmm$A) <- state.labels
	hmm$distributions <- cbind(mean=hmm$means, variance=hmm$variances, size=fsize(hmm$means,hmm$variances), prob=fprob(hmm$means,hmm$variances))
	rownames(hmm$distributions) <- state.labels
	hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A.initial) <- state.labels
	colnames(hmm$A.initial) <- state.labels
	hmm$distributions.initial <- cbind(mean=hmm$means.initial, variance=hmm$variances.initial, size=fsize(hmm$means.initial,hmm$variances.initial), prob=fprob(hmm$means.initial,hmm$variances.initial))
	rownames(hmm$distributions.initial) <- state.labels

	# Delete redundant entries
	hmm$r <- NULL
	hmm$p <- NULL
	hmm$r.initial <- NULL
	hmm$p.initial <- NULL

	# Issue warnings
	if (num.trials == 1) {
		if (hmm$loglik.delta > hmm$eps) {
			war <- warning("HMM did not converge!\n")
		}
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
