findCNVs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE) {

	## Intercept user input
	IDcheck <- ID  #trigger error if not defined
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
	read.cutoff <- as.integer(quantile(reads, 0.9999))
	if (filter.reads) {
		mask <- reads > read.cutoff
		reads[mask] <- read.cutoff
		numfiltered <- length(which(mask))
		cat(paste0("Replaced read counts > ",read.cutoff," (99.99% quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'filter.reads=FALSE' to disable this filtering.\n"))
# 		if (numfiltered > 0) {
# 			warning(paste("There are very high read counts in your data (probably artificial). Replaced read counts > ",read.cutoff," (99.99% quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'filter.reads=FALSE' to disable this filtering.", sep=""))
# 		}
	}
	
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			size = double(length=numstates), # double* size
			prob = double(length=numstates), # double* prob
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
			size.initial = double(length=numstates), # double* initial_size
			prob.initial = double(length=numstates), # double* initial_prob
			A.initial = double(length=numstates*numstates), # double* initial_A
			proba.initial = double(length=numstates), # double* initial_proba
			use.initial.params = as.logical(0), # bool* use_initial_params
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
			A.initial = as.vector(hmm$A), # double* initial_A
			proba.initial = as.vector(hmm$proba), # double* initial_proba
			use.initial.params = as.logical(1), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)
	}

	# Add useful entries
	hmm$ID <- ID
	names(hmm$weights) <- state.labels
	hmm$coordinates <- data.frame(as.character(seqnames(binned.data)), start(ranges(binned.data)), end(ranges(binned.data)))
	names(hmm$coordinates) <- coordinate.names
	hmm$seqlengths <- seqlengths(binned.data)
	class(hmm) <- class.aneufinder.hmm
	hmm$states <- factor(state.labels, levels=state.labels)[hmm$states+1]
	hmm$eps <- eps
	hmm$A <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A) <- state.labels
	colnames(hmm$A) <- state.labels
	hmm$distributions <- data.frame(type=state.distributions, size=hmm$size, prob=hmm$prob, mu=fmean(hmm$size,hmm$prob), variance=fvariance(hmm$size,hmm$prob))
	rownames(hmm$distributions) <- state.labels
	# Treat 'null-mixed' separately
	if ('null-mixed' %in% state.labels) {
		hmm$distributions['null-mixed','mu'] <- (1-hmm$distributions['null-mixed','prob'])/hmm$distributions['null-mixed','prob']
		hmm$distributions['null-mixed','variance'] <- hmm$distributions['null-mixed','mu']/hmm$distributions['null-mixed','prob']
		hmm$distributions['null-mixed','size'] <- NA
	}
	hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A.initial) <- state.labels
	colnames(hmm$A.initial) <- state.labels
	hmm$distributions.initial <- data.frame(type=state.distributions, size=hmm$size.initial, prob=hmm$prob.initial, mu=fmean(hmm$size.initial,hmm$prob.initial), variance=fvariance(hmm$size.initial,hmm$prob.initial))
	rownames(hmm$distributions.initial) <- state.labels
	hmm$distributions.initial['nullsomy',2:5] <- c(0,1,0,0)
	hmm$filter.reads <- filter.reads

	# Delete redundant entries
	hmm$size <- NULL
	hmm$prob <- NULL
	hmm$size.initial <- NULL
	hmm$prob.initial <- NULL
	hmm$use.initial.params <- NULL
	hmm$read.cutoff <- NULL
	hmm$distr.type <- NULL

	# Issue warnings
	if (num.trials == 1) {
		if (hmm$loglik.delta > hmm$eps) {
			war <- warning("HMM did not converge!\n")
		}
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
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
