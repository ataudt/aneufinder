univariate.from.binned.data = function(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE) {

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


	war = NULL
	if (is.null(eps.try)) eps.try = eps

	coordinatenames = c("chrom","start","end","reads")
	names(binned.data) = coordinatenames

	## Assign variables
	# number of states
	numstates = 3
	# prepare densities to use
	densityNames = c("ZI", "NB", "NB")
	enum = c("ZI","NB","ZINB","Other")
	densityNames = factor(densityNames, levels=enum)
	densityNames = as.numeric(densityNames) - 1
	# length of observation sequence
	numbins = length(binned.data$reads) 

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
	modellist = list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")
		z = .C("R_univariate_hmm",
			as.integer(binned.data$reads), # double* O
			as.integer(numbins), # int* T
			as.integer(numstates), # int* N
			as.integer(densityNames), # I pass them but they are fixed
			double(length=numstates), # double* r
			double(length=numstates), # double* p
			as.integer(max.it), #  int* maxiter
			as.integer(max.time), # double* maxtime
			as.double(eps.try), # double* eps
			double(length=numbins * numstates), # double* post
			double(length=numstates*numstates), # double* A
			double(length=numstates), # double* proba
			double(length=1), # double* loglik
			double(length=numstates), # double* softweights
			double(length=numstates), # double* initial r
			double(length=numstates), # double* initial p
			double(length=numstates*numstates), # double* initial A
			double(length=numstates), # double* initial proba
			as.logical(0), # bool* use_initial_params
			as.integer(num.threads) # int* num_threads
		)

		model = NULL
		model$loglik = z[[13]]
		model$iteration = z[[7]]
		model$time.in.sec = z[[8]]
		model$delta.loglik = z[[9]]
		model$epsilon = eps.try
		model$proba = z[[12]]
		model$A = matrix(z[[11]], ncol=numstates, byrow=TRUE)
		model$distributions = cbind(r=z[[5]], p=z[[6]])
		model$softweights = z[[14]]
		model$proba.initial = z[[18]]
		model$A.initial = matrix(z[[17]], ncol=numstates, byrow=TRUE)
		model$distributions.initial = cbind(r=z[[15]], p=z[[16]])
		model$num.threads = z[[20]]
		if (num.trials > 1) {
			if (model$delta.loglik > model$epsilon) {
				warning("HMM did not converge in trial run ",i_try,"!\n")
			}
			# Store model in list
			modellist[[i_try]] = model
		}
	}

	if (num.trials > 1) {

		# Select fit with best loglikelihood
		indexmax = which.max(unlist(lapply(modellist,"[[","loglik")))
		model = modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		cat("\n\nRerunning try ",indexmax," with eps =",eps,"--------------------\n")
		z = .C("R_univariate_hmm",
			as.integer(binned.data$reads), # double* O
			as.integer(numbins), # int* T
			as.integer(numstates), # int* N
			as.integer(densityNames), # I pass them but they are fixed
			double(length=numstates), # double* r
			double(length=numstates), # double* p
			as.integer(max.it), #  int* maxiter
			as.integer(max.time), # double* maxtime
			as.double(eps), # double* eps
			double(length=numbins * numstates), # double* post
			double(length=numstates*numstates), # double* A
			double(length=numstates), # double* proba
			double(length=1), # double* loglik
			double(length=numstates), # double* softweights
			as.vector(model$distributions[,'r']), # double* initial r
			as.vector(model$distributions[,'p']), # double* initial p
			as.vector(model$A), # double* initial A
			as.vector(model$proba), # double* initial proba
			as.logical(1), # bool* use_initial_params
			as.integer(num.threads) # int* num_threads
		)
	}

	model = NULL
	model$coordinates = binned.data[,coordinatenames[1:3]]
	model$reads = binned.data$reads
	model$posteriors = matrix(z[[10]], ncol=numstates)
	colnames(model$posteriors) = c("P(state = zero inflation)","P(state = unmodified)","P(state = modified)")
	class(model) = "chromStar.univariate.model"
	model$states = get.states(model)
	model$states.with.zeroinflation = get.states(model, separate.zeroinflation=TRUE)
	model$loglik = z[[13]]
	model$iteration = z[[7]]
	model$time.in.sec = z[[8]]
	model$delta.loglik = z[[9]]
	model$epsilon = eps
	model$proba = z[[12]]
	model$A = matrix(z[[11]], ncol=numstates, byrow=TRUE)
	model$distributions = cbind(r=z[[5]], p=z[[6]])
	model$softweights = z[[14]]
	model$proba.initial = z[[18]]
	model$A.initial = matrix(z[[17]], ncol=numstates, byrow=TRUE)
	model$distributions.initial = cbind(r=z[[15]], p=z[[16]])
	model$num.threads = z[[20]]
	if (num.trials == 1) {
		if (model$delta.loglik > model$epsilon) {
			war = warning("HMM did not converge!\n")
		}
	}

	if (!is.null(war)) {
		if (output.if.not.converged == TRUE) {
			return(model)
		}
	} else {
		return(model)
	}
}
