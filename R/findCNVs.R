

#' Find copy number variations
#'
#' \code{findCNVs} classifies the binned read counts into several states which represent copy-number-variation.
#'
#' \code{findCNVs} uses a 6-state Hidden Markov Model to classify the binned read counts: state '0-somy' with a delta function as emission densitiy (only zero read counts), '1-somy','2-somy','3-somy','4-somy', etc. with negative binomials (see \code{\link{dnbinom}}) as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper \code{citation("AneuFinder")} for a detailed description of the method.
#' @author Aaron Taudt
#' @inheritParams univariate.findCNVs
#' @param method One of \code{c('univariate','bivariate','DNAcopy')}. In the univariate case strand information is discarded, while in the bivariate case strand information is used for the fitting. For \code{method="dnacopy"} the \pkg{\link{DNAcopy}} package is used which gives more robust but less sensitive results.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom stats dgeom dnbinom
#' @export
#'
#'@examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the BAM file into bin size 1Mp
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'## Fit the Hidden Markov Model
#'model <- findCNVs(binned[[1]], eps=0.1, max.time=60)
#'## Check the fit
#'plot(model, type='histogram')
#'
findCNVs <- function(binned.data, ID=NULL, eps=0.1, init="standard", max.time=-1, max.iter=1000, num.trials=15, eps.try=10*eps, num.threads=1, count.cutoff.quantile=0.999, strand='*', states=c("zero-inflation",paste0(0:10,"-somy")), most.frequent.state="2-somy", method="univariate", algorithm="EM", initial.params=NULL) {

	## Intercept user input
  binned.data <- loadFromFiles(binned.data, check.class='GRanges')
	if (is.null(ID)) {
		ID <- attr(binned.data, 'ID')
	}

	## Print some stuff
	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	message("Find CNVs for ID = ",ID, ":")

	if (method == 'univariate') {
		model <- univariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, count.cutoff.quantile=count.cutoff.quantile, strand=strand, states=states, most.frequent.state=most.frequent.state, algorithm=algorithm, initial.params=initial.params)
	} else if (method == 'bivariate') {
		model <- bivariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, count.cutoff.quantile=count.cutoff.quantile, states=states, most.frequent.state=most.frequent.state, initial.params=initial.params)
	} else if (method == 'dnacopy') {
	  model <- DNAcopy.findCNVs(binned.data, ID, most.frequent.state='2-somy', count.cutoff.quantile=count.cutoff.quantile, strand=strand)
	}

	attr(model, 'call') <- call
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(model)

}


#' Find copy number variations (univariate)
#'
#' \code{findCNVs} classifies the binned read counts into several states which represent copy-number-variation.
#'
#' @param binned.data A \link{GRanges} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various downstream functions. Could be the file name of the \code{binned.data} for example.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param init One of the following initialization procedures:
#'	\describe{
#'		\item{\code{standard}}{The negative binomial of state '2-somy' will be initialized with \code{mean=mean(counts)}, \code{var=var(counts)}. This procedure usually gives good convergence.}
#'		\item{\code{random}}{Mean and variance of the negative binomial of state '2-somy' will be initialized with random values (in certain boundaries, see source code). Try this if the \code{standard} procedure fails to produce a good fit.}
#'	}
#' @param max.time The maximum running time in seconds for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.
#' @param max.iter The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.
#' @param num.trials The number of trials to find a fit where state \code{most.frequent.state} is most frequent. Each time, the HMM is seeded with different random initial values.
#' @param eps.try If code num.trials is set to greater than 1, \code{eps.try} is used for the trial runs. If unset, \code{eps} is used.
#' @param num.threads Number of threads to use. Setting this to >1 may give increased performance.
#' @param count.cutoff.quantile A quantile between 0 and 1. Should be near 1. Read counts above this quantile will be set to the read count specified by this quantile. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. Set \code{count.cutoff.quantile=1} in this case.
#' @param strand Run the HMM only for the specified strand. One of \code{c('+', '-', '*')}.
#' @param states A subset or all of \code{c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy",...)}. This vector defines the states that are used in the Hidden Markov Model. The order of the entries must not be changed.
#' @param most.frequent.state One of the states that were given in \code{states} or \code{NULL}. The specified state is assumed to be the most frequent one. This can help the fitting procedure to converge into the correct fit.
#' @param algorithm One of \code{c('baumWelch','EM')}. The expectation maximization (\code{'EM'}) will find the most likely states and fit the best parameters to the data, the \code{'baumWelch'} will find the most likely states using the initial parameters.
#' @param initial.params A \code{\link{aneuHMM}} object or file containing such an object from which initial starting parameters will be extracted.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom stats runif
univariate.findCNVs <- function(binned.data, ID=NULL, eps=0.1, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, count.cutoff.quantile=0.999, strand='*', states=c("zero-inflation",paste0(0:10,"-somy")), most.frequent.state="2-somy", algorithm="EM", initial.params=NULL) {

	### Define cleanup behaviour ###
	on.exit(.C("C_univariate_cleanup", PACKAGE = 'AneuFinder'))

	## Intercept user input
  binned.data <- loadFromFiles(binned.data, check.class='GRanges')[[1]]
	if (is.null(ID)) {
		ID <- attr(binned.data, 'ID')
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.strand(strand)!=0) stop("argument 'strand' expects either '+', '-' or '*'")
  if (!is.null(most.frequent.state)) {
    	if (!most.frequent.state %in% states) stop("argument 'most.frequent.state' must be one of c(",paste(states, collapse=","),") or NULL")
  }
	if (!algorithm %in% c('baumWelch','EM')) {
		stop("argument 'algorithm' expects one of c('baumWelch','EM')")
	}
	if (algorithm == 'baumWelch' & num.trials>1) {
		warning("Set 'num.trials <- 1' because 'algorithm==\"baumWelch\"'.")
		num.trials <- 1
	}
	initial.params <- loadFromFiles(initial.params, check.class=class.univariate.hmm)[[1]]
	if (class(initial.params)!=class.univariate.hmm & !is.null(initial.params)) {
		stop("argument 'initial.params' expects a ",class.univariate.hmm," object or file that contains such an object")
	}
	if (algorithm == 'baumWelch' & is.null(initial.params)) {
		warning("'initial.params' should be specified if 'algorithm=\"baumWelch\"")
	}
	if (!is.null(initial.params)) {
		init <- 'initial.params'
	}

	warlist <- list()
	if (num.trials==1) eps.try <- eps

	## Assign variables
	inistates <- initializeStates(states)
	state.labels <- inistates$states
	state.distributions <- inistates$distributions
	multiplicity <- inistates$multiplicity
	dependent.states.mask <- (state.labels != 'zero-inflation') & (state.labels != '0-somy')
	numstates <- length(states)
	numbins <- length(binned.data)
	iniproc <- which(init==c("standard","random")) # transform to int
	if (strand=='+') {
		select <- 'pcounts'
	} else if (strand=='-') {
		select <- 'mcounts'
	} else if (strand=='*') {
		select <- 'counts'
	}
	counts <- mcols(binned.data)[,select]
	algorithm <- factor(algorithm, levels=c('baumWelch','viterbi','EM'))

	### Make return object
		result <- list()
		class(result) <- class.univariate.hmm
		result$ID <- ID
		result$bins <- binned.data
	## Quality info
		qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=NA)
		result$qualityInfo <- qualityInfo
	## Convergence info
		convergenceInfo <- list(eps=eps, loglik=NA, loglik.delta=NA, num.iterations=NA, time.sec=NA, error=NA)
		result$convergenceInfo <- convergenceInfo

	# Check if there are counts in the data, otherwise HMM will blow up
	if (any(is.na(counts))) {
		stop(paste0("ID = ",ID,": NAs found in reads."))
	}
	if (!any(counts!=0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No HMM done."))
		result$warnings <- warlist
		return(result)
	} else if (any(counts<0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No HMM done."))
		result$warnings <- warlist
		return(result)
	}
		

	# Filter high counts out, makes HMM faster
	count.cutoff <- quantile(counts, count.cutoff.quantile)
	names.count.cutoff <- names(count.cutoff)
	count.cutoff <- ceiling(count.cutoff)
	mask <- counts > count.cutoff
	counts[mask] <- count.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		message(paste0("Replaced read counts > ",count.cutoff," (",names.count.cutoff," quantile) by ",count.cutoff," in ",numfiltered," bins. Set option 'count.cutoff.quantile=1' to disable this filtering. This filtering was done to enhance performance."))
	}
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		message(paste0("Trial ",i_try," / ",num.trials))

		## Initial parameters
		if (init == 'initial.params') {
			A.initial <- initial.params$transitionProbs
			proba.initial <- initial.params$startProbs
			size.initial <- initial.params$distributions[,'size']
			prob.initial <- initial.params$distributions[,'prob']
			size.initial[is.na(size.initial)] <- 0
			prob.initial[is.na(prob.initial)] <- 0
		} else if (init == 'random') {
			A.initial <- matrix(stats::runif(numstates^2), ncol=numstates)
			A.initial <- sweep(A.initial, 1, rowSums(A.initial), "/")			
			proba.initial <- stats::runif(numstates)
			# Distributions for dependent states
			size.initial <- stats::runif(1, min=0, max=100) * cumsum(dependent.states.mask)
			prob.initial <- stats::runif(1) * dependent.states.mask
			# Assign initials for the 0-somy distribution
			index <- which('0-somy'==state.labels)
			size.initial[index] <- 1
			prob.initial[index] <- 0.5
		} else if (init == 'standard') {
			A.initial <- matrix(NA, ncol=numstates, nrow=numstates)
			for (irow in 1:numstates) {
				for (icol in 1:numstates) {
					if (irow==icol) { A.initial[irow,icol] <- 0.9 }
					else { A.initial[irow,icol] <- 0.1/(numstates-1) }
				}
			}
			proba.initial <- rep(1/numstates, numstates)
			## Set initial mean of most.frequent.state distribution to max of count histogram
			max.counts <- as.integer(names(which.max(table(counts[counts>0]))))
			divf <- max(multiplicity[most.frequent.state], 1)
			mean.initial.monosomy <- max.counts/divf
			var.initial.monosomy <- var(counts) / 2
# 			mean.initial.monosomy <- mean(counts[counts>0])/divf
# 			var.initial.monosomy <- var(counts[counts>0])/divf
			if (is.na(mean.initial.monosomy)) {
				mean.initial.monosomy <- 1
			}
			if (is.na(var.initial.monosomy)) {
				var.initial.monosomy <- mean.initial.monosomy + 1
			}
			if (mean.initial.monosomy >= var.initial.monosomy) {
				mean.initial <- mean.initial.monosomy * cumsum(dependent.states.mask)
				var.initial <- (mean.initial.monosomy+1) * cumsum(dependent.states.mask)
				size.initial <- rep(0,numstates)
				prob.initial <- rep(0,numstates)
				mask <- dependent.states.mask
				size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
				prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
			} else {
				mean.initial <- mean.initial.monosomy * cumsum(dependent.states.mask)
				var.initial <- var.initial.monosomy * cumsum(dependent.states.mask)
				size.initial <- rep(0,numstates)
				prob.initial <- rep(0,numstates)
				mask <- dependent.states.mask
				size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
				prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
			}
			# Assign initials for the 0-somy distribution
			index <- which('0-somy'==state.labels)
			size.initial[index] <- 1
			prob.initial[index] <- 0.5
		}
	
		hmm <- .C("C_univariate_hmm",
			counts = as.integer(counts), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			state.labels = as.integer(state.labels), # int* state_labels
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
			size.initial = as.vector(size.initial), # double* initial_size
			prob.initial = as.vector(prob.initial), # double* initial_prob
			A.initial = as.vector(A.initial), # double* initial_A
			proba.initial = as.vector(proba.initial), # double* initial_proba
			use.initial.params = as.logical(1), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			count.cutoff = as.integer(count.cutoff), # int* count.cutoff
			algorithm = as.integer(algorithm), # int* algorithm
			PACKAGE = 'AneuFinder'
		)

		hmm$eps <- eps.try
		if (num.trials > 1) {
			if (hmm$loglik.delta > hmm$eps) {
				warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": HMM did not converge in trial run ",i_try,"!\n"))
			}
			# Store model in list
			hmm$counts <- NULL
			modellist[[as.character(i_try)]] <- hmm
			init <- 'random'
		} else if (num.trials == 1) {
			if (hmm$loglik.delta > eps) {
				warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": HMM did not converge!\n"))
			}
		}
	}

	if (num.trials > 1) {

		# Mathematically we should select the fit with highest loglikelihood. If we think the fit with the highest loglikelihood is incorrect, we should change the underlying model. However, this is very complex and we choose to select a fit that we think is (more) correct, although it has not the highest support given our (imperfect) model.
		if (length(modellist)>1) {
			## Select models where weight of most.frequent.state is at least half of that of actual most frequent state, then select model with highest loglik
			logliks <- unlist(lapply(modellist,'[[','loglik'))
			df.weight <- as.data.frame(lapply(modellist, '[[', 'weights'))
			names(df.weight) <- 1:length(modellist)
			rownames(df.weight) <- state.labels
			models2use <- df.weight[most.frequent.state,] / apply(df.weight, 2, max) > 0.5
			models2use[is.na(models2use)] <- FALSE
			if (any(models2use)) {
				index2use <- names(which.max(logliks[models2use]))
			} else {
				index2use <- names(which.max(logliks))
			}
		} else {
			index2use <- 1
		}
		hmm <- modellist[[index2use]]

		# Check if size and prob parameter are correct
		if (any(is.na(hmm$size) | is.nan(hmm$size) | is.infinite(hmm$size) | is.na(hmm$prob) | is.nan(hmm$prob) | is.infinite(hmm$prob))) {
			hmm$error <- 3
		} else {

			# Rerun the HMM with different epsilon and initial parameters from trial run
			message(paste0("Rerunning trial ",index2use," with eps = ",eps))
			hmm <- .C("C_univariate_hmm",
				counts = as.integer(counts), # int* O
				num.bins = as.integer(numbins), # int* T
				num.states = as.integer(numstates), # int* N
				state.labels = as.integer(state.labels), # int* state_labels
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
				size.initial = as.vector(hmm$size), # double* initial_size
				prob.initial = as.vector(hmm$prob), # double* initial_prob
				A.initial = as.vector(hmm$A), # double* initial_A
				proba.initial = as.vector(hmm$proba), # double* initial_proba
				use.initial.params = as.logical(1), # bool* use_initial_params
				num.threads = as.integer(num.threads), # int* num_threads
				error = as.integer(0), # int* error (error handling)
				count.cutoff = as.integer(count.cutoff), # int* count.cutoff
				algorithm = as.integer(algorithm), # int* algorithm
				PACKAGE = 'AneuFinder'
			)
		}

	}

	### Make return object ###
	## Check for errors
		if (hmm$error == 0) {
		## Bin coordinates and states ###
			result$bins$state <- state.labels[hmm$states]
			result$bins$copy.number <- multiplicity[hmm$states]
		## Segmentation
			message("Making segmentation ...", appendLF=FALSE)
			ptm <- proc.time()
			suppressMessages(
				result$segments <- as(collapseBins(as.data.frame(result$bins), column2collapseBy='state', columns2drop='width', columns2average=c('counts','mcounts','pcounts')), 'GRanges')
			)
			seqlevels(result$segments) <- seqlevels(result$bins) # correct order from as()
			seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
			time <- proc.time() - ptm
			message(" ",round(time[3],2),"s")
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
			names(result$startProbs) <- state.labels
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- state.labels
			# Distributions
				distributions <- data.frame()
				distributions.initial <- data.frame()
				for (idistr in 1:length(hmm$distr.type)) {
					distr <- levels(state.distributions)[hmm$distr.type[idistr]]
					if (distr == 'dnbinom') {
						distributions <- rbind(distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], mu=dnbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dnbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
						distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], mu=dnbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dnbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
					} else if (distr == 'dgeom') {
						distributions <- rbind(distributions, data.frame(type=distr, size=NA, prob=hmm$prob[idistr], mu=dgeom.mean(hmm$prob[idistr]), variance=dgeom.variance(hmm$prob[idistr])))
						distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=NA, prob=hmm$prob.initial[idistr], mu=dgeom.mean(hmm$prob.initial[idistr]), variance=dgeom.variance(hmm$prob.initial[idistr])))
					} else if (distr == 'delta') {
						distributions <- rbind(distributions, data.frame(type=distr, size=NA, prob=NA, mu=0, variance=0))
						distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=NA, prob=NA, mu=0, variance=0))
					} else if (distr == 'dbinom') {
						distributions <- rbind(distributions, data.frame(type=distr, size=hmm$size[idistr], prob=hmm$prob[idistr], mu=dbinom.mean(hmm$size[idistr],hmm$prob[idistr]), variance=dbinom.variance(hmm$size[idistr],hmm$prob[idistr])))
						distributions.initial <- rbind(distributions.initial, data.frame(type=distr, size=hmm$size.initial[idistr], prob=hmm$prob.initial[idistr], mu=dbinom.mean(hmm$size.initial[idistr],hmm$prob.initial[idistr]), variance=dbinom.variance(hmm$size.initial[idistr],hmm$prob.initial[idistr])))
					}
				}
				rownames(distributions) <- state.labels
				rownames(distributions.initial) <- state.labels
				result$distributions <- distributions
				result$distributions.initial <- distributions.initial
		## Convergence info
			convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec, error=hmm$error)
			result$convergenceInfo <- convergenceInfo
		## Quality info
			qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=qc.bhattacharyya(result))
			result$qualityInfo <- qualityInfo
		} else if (hmm$error == 1) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'count.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
		} else if (hmm$error == 2) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": An error occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library"))
		} else if (hmm$error == 3) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": NA/NaN/Inf in 'size' or 'prob' parameter detected. This is probably because your binned data contains too few bins."))
		}

	## Issue warnings
	result$warnings <- warlist

	## Return results
	return(result)
}


#' Find copy number variations (bivariate)
#'
#' \code{bivariate.findCNVs} finds CNVs using read count information from both strands.
#'
#' @inheritParams univariate.findCNVs
#' @return An \code{\link{aneuBiHMM}} object.
#' @importFrom stats pgeom pnbinom qnorm
bivariate.findCNVs <- function(binned.data, ID=NULL, eps=0.1, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, count.cutoff.quantile=0.999, states=c("zero-inflation",paste0(0:10,"-somy")), most.frequent.state="1-somy", algorithm='EM', initial.params=NULL) {

	## Intercept user input
  binned.data <- loadFromFiles(binned.data, check.class='GRanges')[[1]]
	if (is.null(ID)) {
		ID <- attr(binned.data, 'ID')
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	initial.params <- loadFromFiles(initial.params, check.class=class.bivariate.hmm)[[1]]
	if (class(initial.params)!=class.bivariate.hmm & !is.null(initial.params)) {
		stop("argument 'initial.params' expects a ",class.bivariate.hmm," object or file that contains such an object")
	}
	if (algorithm == 'baumWelch' & is.null(initial.params)) {
		warning("'initial.params' should be specified if 'algorithm=\"baumWelch\"")
	}
	if (!is.null(initial.params)) {
		init <- 'initial.params'
	}

	warlist <- list()
	if (is.null(eps.try)) eps.try <- eps

	## Variables
	num.bins <- length(binned.data)
	inistates <-  initializeStates(states)
	multiplicity <- inistates$multiplicity
	state.labels <- inistates$states
	algorithm <- factor(algorithm, levels=c('baumWelch','viterbi','EM'))

	## Get counts
	select <- 'counts'
	counts <- matrix(c(mcols(binned.data)[,paste0('m',select)], mcols(binned.data)[,paste0('p',select)]), ncol=2, dimnames=list(bin=1:num.bins, strand=c('minus','plus')))
	maxcounts <- max(counts)

	## Filter high counts out, makes HMM faster
	count.cutoff <- quantile(counts, count.cutoff.quantile)
	names.count.cutoff <- names(count.cutoff)
	count.cutoff <- ceiling(count.cutoff)
	mask <- counts > count.cutoff
	counts[mask] <- count.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		message(paste0("Replaced read counts > ",count.cutoff," (",names.count.cutoff," quantile) by ",count.cutoff," in ",numfiltered," bins. Set option 'count.cutoff.quantile=1' to disable this filtering. This filtering was done to enhance performance."))
	}
	
	### Make return object
		result <- list()
		class(result) <- class.bivariate.hmm
		result$ID <- ID
		result$bins <- binned.data
	## Quality info
		qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=NA)
		result$qualityInfo <- qualityInfo

	# Check if there are counts in the data, otherwise HMM will blow up
	if (!any(counts!=0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No HMM done."))
		result$warnings <- warlist
		return(result)
	} else if (any(counts<0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No HMM done."))
		result$warnings <- warlist
		return(result)
	}

	if (init=='initial.params') {
		uni.transitionProbs <- initial.params$univariateParams$transitionProbs
		uni.startProbs <- initial.params$univariateParams$startProbs
		distributions <- initial.params$distributions
		uni.weights <- initial.params$univariateParams$weights
		uni.states <- names(uni.weights)
		num.uni.states <- length(uni.states)
		num.models <- length(distributions)
		num.bins <- length(initial.params$bins)
		comb.states <- factor(names(initial.params$startProbs), levels=names(initial.params$startProbs))
		num.comb.states <- length(comb.states)
		states.list <- list(minus=initial.params$bins$mstate, plus=initial.params$bins$pstate)
		comb.states.per.bin <- factor(do.call(paste, states.list), levels=levels(comb.states))
		A.initial <- initial.params$transitionProbs
		proba.initial <- initial.params$startProbs
		use.initial <- TRUE
	} else {
		### Stack the strands and run one univariate findCNVs
		message("")
		message(paste(rep('-',getOption('width')), collapse=''))
		message("Running univariate")
		binned.data.minus <- binned.data
		strand(binned.data.minus) <- '-'
		binned.data.minus$counts <- binned.data.minus$mcounts
		binned.data.plus <- binned.data
		strand(binned.data.plus) <- '+'
		binned.data.plus$counts <- binned.data.plus$pcounts
		binned.data.stacked <- c(binned.data.minus, binned.data.plus)
		mask.attributes <- c('qualityInfo', 'ID', 'min.mapq')
		attributes(binned.data.stacked)[mask.attributes] <- attributes(binned.data)[mask.attributes]

		model.stacked <- univariate.findCNVs(binned.data.stacked, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, count.cutoff.quantile=1, states=states, most.frequent.state=most.frequent.state)
		model.minus <- model.stacked
		model.minus$bins <- model.minus$bins[strand(model.minus$bins)=='-']
		model.minus$segments <- model.minus$segments[strand(model.minus$segments)=='-']
		model.plus <- model.stacked
		model.plus$bins <- model.plus$bins[strand(model.plus$bins)=='+']
		model.plus$segments <- model.plus$segments[strand(model.plus$segments)=='+']

		models <- list(minus=model.minus, plus=model.plus)

		## Extract counts and other stuff
		uni.transitionProbs <- lapply(models, '[[', 'transitionProbs')[[1]]
		uni.startProbs <- lapply(models, '[[', 'startProbs')[[1]]
		distributions <- lapply(models, '[[', 'distributions')
		uni.weights <- lapply(models, '[[', 'weights')[[1]]
		uni.states <- names(uni.weights)
		num.uni.states <- length(uni.states)
		num.models <- length(distributions)
		num.bins <- length(models[[1]]$bins)
		comb.states <- vector()
		levels.state <- levels(models[[1]]$bins$state)
		for (i1 in 1:length(levels.state)) {
			for (i2 in 1:length(levels.state)) {
				comb.state <- paste(levels.state[i1], levels.state[i2])
				comb.states[length(comb.states)+1] <- comb.state
			}
		}
		comb.states <- factor(comb.states, levels=comb.states)
		num.comb.states <- length(comb.states)
		states.list <- lapply(models, function(x) { x$bins$state })
		comb.states.per.bin <- factor(do.call(paste, states.list), levels=levels(comb.states))
		A.initial <- double(length=num.comb.states*num.comb.states)
		proba.initial <- double(length=num.comb.states)
		use.initial <- FALSE
	}

	### Prepare the multivariate HMM
	message("")
	message(paste(rep('-',getOption('width')), collapse=''))
	message("Preparing bivariate HMM\n")

	## Pre-compute z-values for each number of counts
	ptm <- startTimedMessage("Computing pre z-matrix...")
	z.per.count <- array(NA, dim=c(maxcounts+1, num.models, num.uni.states), dimnames=list(counts=0:maxcounts, strand=names(distributions), state=uni.states))
	xcounts <- 0:maxcounts
	for (istrand in 1:num.models) {
		for (istate in 1:num.uni.states) {
			if (distributions[[istrand]][istate,'type']=='dnbinom') {
				size <- distributions[[istrand]][istate,'size']
				prob <- distributions[[istrand]][istate,'prob']
				u <- stats::pnbinom(xcounts, size, prob)
			} else if (distributions[[istrand]][istate,'type']=='delta') {
				u <- rep(1, length(xcounts))
			} else if (distributions[[istrand]][istate,'type']=='dgeom') {
				prob <- distributions[[istrand]][istate,'prob']
				u <- stats::pgeom(xcounts, prob)
			}
			qnorm_u <- stats::qnorm(u)
			mask <- qnorm_u==Inf
			qnorm_u[mask] <- stats::qnorm(1-1e-16)
			z.per.count[ , istrand, istate] <- qnorm_u
		}
	}
	stopTimedMessage(ptm)

	## Compute the z matrix
	ptm <- startTimedMessage("Transfering values into z-matrix...")
	z.per.bin = array(NA, dim=c(num.bins, num.models, num.uni.states), dimnames=list(bin=1:num.bins, strand=names(distributions), state=uni.states))
	for (istrand in 1:num.models) {
		for (istate in 1:num.uni.states) {
			z.per.bin[ , istrand, istate] = z.per.count[counts[,istrand]+1, istrand, istate]
		}
	}
	remove(z.per.count)
	stopTimedMessage(ptm)

	## Calculate correlation matrix
	ptm <- startTimedMessage("Computing inverse of correlation matrix...")
	correlationMatrix <- array(0, dim=c(num.models,num.models,num.comb.states), dimnames=list(strand=names(distributions), strand=names(distributions), comb.state=comb.states))
	correlationMatrixInverse <- correlationMatrix
	determinant <- rep(0, num.comb.states)
	names(determinant) <- comb.states
	usestateTF <- rep(NA, num.comb.states) # TRUE, FALSE vector for usable states
	names(usestateTF) <- comb.states
	for (comb.state in comb.states) {
		state <- strsplit(comb.state, ' ')[[1]]
		mask <- which(comb.states.per.bin==comb.state)
		# Subselect z
		z.temp <- matrix(NA, ncol=num.models, nrow=length(mask))
		for (istrand in 1:num.models) {
			z.temp[,istrand] <- z.per.bin[mask, istrand, state[istrand]]
		}
		temp <- tryCatch({
			if (nrow(z.temp) > 1) {
				correlationMatrix[,,comb.state] <- cor(z.temp)
				determinant[comb.state] <- det( correlationMatrix[,,comb.state] )
				correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
				usestateTF[comb.state] <- TRUE
			} else {
				correlationMatrix[,,comb.state] <- diag(num.models)
				determinant[comb.state] <- det( correlationMatrix[,,comb.state] )
				correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
				usestateTF[comb.state] <- TRUE
			}
			0
		}, warning = function(war) {
			1
		}, error = function(err) {
			1
		})
		if (temp!=0) {
			correlationMatrix[,,comb.state] <- diag(num.models)
			determinant[comb.state] <- det( correlationMatrix[,,comb.state] )
			correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
			usestateTF[comb.state] <- TRUE
		}
	}
	stopTimedMessage(ptm)

	# Use only states with valid correlationMatrixInverse (all states are valid in this version)
	correlationMatrix = correlationMatrix[,,usestateTF]
	correlationMatrixInverse = correlationMatrixInverse[,,usestateTF]
	comb.states = comb.states[usestateTF]
	comb.states <- droplevels(comb.states)
	determinant = determinant[usestateTF]
	num.comb.states <- length(comb.states)

	### Calculate multivariate densities for each state ###
	ptm <- startTimedMessage("Calculating multivariate densities...")
	densities <- matrix(1, ncol=num.comb.states, nrow=num.bins, dimnames=list(bin=1:num.bins, comb.state=comb.states))
	for (comb.state in comb.states) {
		istate <- which(comb.state==comb.states)
		state <- strsplit(comb.state, ' ')[[1]]
		z.temp <- matrix(NA, ncol=num.models, nrow=num.bins)
		product <- 1
		for (istrand in 1:num.models) {
			z.temp[,istrand] <- z.per.bin[, istrand, state[istrand]]
			if (distributions[[istrand]][state[istrand],'type'] == 'dnbinom') {
				size <- distributions[[istrand]][state[istrand],'size']
				prob <- distributions[[istrand]][state[istrand],'prob']
				product <- product * stats::dnbinom(counts[,istrand], size, prob)
			} else if (distributions[[istrand]][state[istrand],'type'] == 'dgeom') {
				prob <- distributions[[istrand]][state[istrand],'prob']
				product <- product * stats::dgeom(counts[,istrand], prob)
			} else if (distributions[[istrand]][state[istrand],'type'] == 'delta') {
				product <- product * ifelse(counts[,istrand]==0, 1, 0)
			}
		}
		exponent <- -0.5 * apply( ( z.temp %*% (correlationMatrixInverse[ , , istate] - diag(num.models)) ) * z.temp, 1, sum)
		densities[,istate] <- product * determinant[istate]^(-0.5) * exp( exponent )
	}
	# Check if densities are > 1
# 	if (any(densities>1)) stop("Densities > 1")
# 	if (any(densities<0)) stop("Densities < 0")
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
	stopTimedMessage(ptm)
		
	### Define cleanup behaviour ###
	on.exit(.C("C_multivariate_cleanup", as.integer(num.comb.states), PACKAGE = 'AneuFinder'))

	### Run the multivariate HMM
	# Call the C function
	hmm <- .C("C_multivariate_hmm",
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
		A.initial = as.vector(A.initial), # double* initial_A
		proba.initial = as.vector(proba.initial), # double* initial_proba
		use.initial.params = as.logical(use.initial), # bool* use_initial_params
		num.threads = as.integer(num.threads), # int* num_threads
		error = as.integer(0), # error handling
		algorithm = as.integer(algorithm), # int* algorithm
		PACKAGE = 'AneuFinder'
		)
			
	### Check convergence ###
	war <- NULL
	if (hmm$loglik.delta > eps) {
		war <- warning(paste0("ID = ",ID,": HMM did not converge!\n"))
	}
	if (hmm$error == 1) {
		warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'count.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
	} else if (hmm$error == 2) {
		warning(paste0("ID = ",ID,": An error occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library"))
	}

	### Make return object ###
	if (hmm$error == 0) {
		result <- list()
		class(result) <- class.bivariate.hmm
		result$ID <- ID
	## Bin coordinates and states
		result$bins <- binned.data
		result$bins$state <- comb.states[hmm$states]
		# Get states as factors in data.frame
		matrix.states <- matrix(unlist(strsplit(as.character(result$bins$state), split=' ')), byrow=TRUE, ncol=num.models, dimnames=list(bin=1:num.bins, strand=names(distributions)))
		names <- c('mstate','pstate')
		for (i1 in 1:num.models) {
			mcols(result$bins)[names[i1]] <- factor(matrix.states[,i1], levels=uni.states)
		}
	## Segmentation
		ptm <- startTimedMessage("Making segmentation ...")
		suppressMessages(
			result$segments <- as(collapseBins(as.data.frame(result$bins), column2collapseBy='state', columns2drop='width', columns2average=c('counts','mcounts','pcounts')), 'GRanges')
		)
		seqlevels(result$segments) <- seqlevels(result$bins) # correct order from as()
		seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]
		stopTimedMessage(ptm)
	## CNV state for both strands combined
		# Bins
		state <- multiplicity[result$bins$mstate] + multiplicity[result$bins$pstate]
		state[state>max(multiplicity)] <- max(multiplicity)
		multiplicity.inverse <- names(multiplicity)
		names(multiplicity.inverse) <- multiplicity
		state <- multiplicity.inverse[as.character(state)]
		state[(result$bins$mstate=='0-somy' | result$bins$pstate=='0-somy') & state=='zero-inflation'] <- '0-somy'
    result$bins$state <- factor(state, levels=names(multiplicity))
		# Segments
		str <- strsplit(as.character(result$segments$state),' ')
		result$segments$mstate <- factor(unlist(lapply(str, '[[', 1)), levels=state.labels)
		result$segments$pstate <- factor(unlist(lapply(str, '[[', 2)), levels=state.labels)
		state <- multiplicity[result$segments$mstate] + multiplicity[result$segments$pstate]
		state[state>max(multiplicity)] <- max(multiplicity)
		multiplicity.inverse <- names(multiplicity)
		names(multiplicity.inverse) <- multiplicity
		state <- multiplicity.inverse[as.character(state)]
		state[(result$segments$mstate=='0-somy' | result$segments$pstate=='0-somy') & state=='zero-inflation'] <- '0-somy'
    result$segments$state <- factor(state, levels=names(multiplicity))
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
			names(result$startProbs) <- comb.states
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- comb.states
			# Distributions
			result$distributions <- distributions
		## Convergence info
			convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec)
			result$convergenceInfo <- convergenceInfo
		## Quality info
			qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=qc.bhattacharyya(result))
			result$qualityInfo <- qualityInfo
		## Univariate infos
			univariateParams <- list(transitionProbs=uni.transitionProbs, startProbs=uni.startProbs, distributions=distributions[[1]], weights=uni.weights)
			result$univariateParams <- univariateParams
	} else if (hmm$error == 1) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'count.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
	} else if (hmm$error == 2) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": An error occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library"))
	}

	## Issue warnings
	result$warnings <- warlist

	## Return results
	return(result)
}


#' Find copy number variations (DNAcopy)
#'
#' \code{findCNVs} classifies the binned read counts into several states which represent copy-number-variation.
#'
#' @param binned.data A \link{GRanges} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various downstream functions. Could be the file name of the \code{binned.data} for example.
#' @param most.frequent.state One of the states that were given in \code{states} or \code{NULL}. The specified state is assumed to be the most frequent one. This can help the fitting procedure to converge into the correct fit.
#' @param count.cutoff.quantile A quantile between 0 and 1. Should be near 1. Read counts above this quantile will be set to the read count specified by this quantile. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. Set \code{count.cutoff.quantile=1} in this case.
#' @param strand Run the HMM only for the specified strand. One of \code{c('+', '-', '*')}.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom DNAcopy CNA smooth.CNA segment
DNAcopy.findCNVs <- function(binned.data, ID=NULL, most.frequent.state='2-somy', count.cutoff.quantile=0.999, strand='*') {

    ## Function definitions
    mean0 <- function(x) {
        y <- x[x>0]
        if (length(y) > 0) {
            return(mean(y))
        } else {
            return(mean(x))
        }
    }
  	## Intercept user input
    binned.data <- loadFromFiles(binned.data, check.class='GRanges')[[1]]
  	if (is.null(ID)) {
    		ID <- attr(binned.data, 'ID')
  	}
  	if (check.strand(strand)!=0) stop("argument 'strand' expects either '+', '-' or '*'")
  
  	warlist <- list()
  
  	## Assign variables
  	if (strand=='+') {
    		select <- 'pcounts'
  	} else if (strand=='-') {
    		select <- 'mcounts'
  	} else if (strand=='*') {
    		select <- 'counts'
  	}
  	counts <- mcols(binned.data)[,select]

  	### Make return object
		result <- list()
		class(result) <- class.univariate.hmm
		result$ID <- ID
		result$bins <- binned.data
  	## Quality info
		qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=NA)
		result$qualityInfo <- qualityInfo

  	# Check if there are counts in the data, otherwise HMM will blow up
  	if (any(is.na(counts))) {
    		stop(paste0("ID = ",ID,": NAs found in reads."))
  	}
  	if (!any(counts!=0)) {
    		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No CNVs found."))
    		result$warnings <- warlist
    		return(result)
  	} else if (any(counts<0)) {
    		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No CNVs found."))
    		result$warnings <- warlist
    		return(result)
  	}
		

  	# Filter high counts out, makes HMM faster
  	count.cutoff <- quantile(counts, count.cutoff.quantile)
  	names.count.cutoff <- names(count.cutoff)
  	count.cutoff <- ceiling(count.cutoff)
  	mask <- counts > count.cutoff
  	counts[mask] <- count.cutoff
  	numfiltered <- length(which(mask))
  	if (numfiltered > 0) {
  		message(paste0("Replaced read counts > ",count.cutoff," (",names.count.cutoff," quantile) by ",count.cutoff," in ",numfiltered," bins. Set option 'count.cutoff.quantile=1' to disable this filtering. This filtering was done to enhance performance."))
  	}
  	
    
  	### DNAcopy ###
  	logcounts <- log2(counts+1 / mean0(counts+1))
    CNA.object <- DNAcopy::CNA(genomdat=logcounts, chrom=as.vector(seqnames(binned.data)), maploc=as.numeric(start(binned.data)), data.type='logratio')
    CNA.smoothed <- DNAcopy::smooth.CNA(CNA.object)
    CNA.segs <- DNAcopy::segment(CNA.smoothed, verbose=0, min.width=5)
    segs <- CNA.segs$output
    segs.splt <- split(segs, segs$chrom)
    for (i1 in 1:length(segs.splt)) {
        segs.splt[[i1]]$loc.end <- c(segs.splt[[i1]]$loc.start[-1] - 1, seqlengths(binned.data)[names(segs.splt)[i1]])
    }
    segs <- do.call(rbind, segs.splt)
  	segs.gr <- GRanges(seqnames=segs$chrom, ranges=IRanges(start=segs$loc.start, end=segs$loc.end))
  	segs.gr$mean.count <- 2^segs$seg.mean
  	
  	## Modify bins to contain median count
  	ind <- findOverlaps(binned.data, segs.gr, select='first')
  	segs.gr$median.count <- sapply(split(counts, ind), median)
    counts.median <- segs.gr$median.count[ind]
    if (mean0(counts.median) > 0) {
        counts.median <- counts.median / mean0(counts.median)
    }
  
    if (is.null(most.frequent.state)) {
        ## Determine Copy Number
        counts.norm <- counts / mean0(counts)
        CNgrid       <- seq(1.5, 6, by=0.05)
        outerRaw     <- counts.norm %o% CNgrid
        outerRound   <- round(outerRaw)
        outerDiff    <- (outerRaw - outerRound) ^ 2
        outerColsums <- colSums(outerDiff, na.rm = FALSE, dims = 1)
        CNmult       <- CNgrid[order(outerColsums)]
        CNerror      <- round(sort(outerColsums), digits=2)
        CN <- CNmult[1]
        # plot(CNgrid, outerColsums)
    } else {
        ## Determine copy numbers such that most.frequent.state is the most frequent
        mfs <- as.integer(sub('-somy','', most.frequent.state))
        CNgrid       <- seq(1.5, 6, by=0.05)
        outerCNstates <- round(outer(counts.median, CNgrid))
        freqs <- apply(outerCNstates, 2, table)
        index <- which.max(sapply(freqs, function(x) { x[as.character(mfs)] }))
        CN <- CNgrid[index]
    }
    
    CN.states <- round(counts.median * CN)
    somies <- paste0(CN.states, '-somy')
    inistates <- suppressWarnings( initializeStates(paste0(sort(unique(CN.states)), '-somy')) )
  	state.labels <- inistates$states
  	state.distributions <- inistates$distributions
  	multiplicity <- inistates$multiplicity
  	
  	
  	### Make return object ###
    result$bins$state <- factor(somies, levels=inistates$states)
  	## Segmentation
		message("Making segmentation ...", appendLF=FALSE)
		ptm <- proc.time()
		suppressMessages(
  			result$segments <- as(collapseBins(as.data.frame(result$bins), column2collapseBy='state', columns2drop='width', columns2average=c('counts','mcounts','pcounts')), 'GRanges')
		)
		seqlevels(result$segments) <- seqlevels(result$bins) # correct order from as()
		seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
  	## Parameters
		# Weights
		tab <- table(result$bins$state)
		result$weights <- tab / sum(tab)
		# Distributions
		distributions <- list()
		bins.splt <- split(result$bins, result$bins$state)
		for (i1 in 1:length(bins.splt)) {
		    qus <- quantile(bins.splt[[i1]]$counts, c(0.01, 0.99))
		    qcounts <- bins.splt[[i1]]$counts
		    qcounts <- qcounts[qcounts >= qus[1] & qcounts <= qus[2]]
		    if (length(qcounts) == 0) {
		        qcounts <- bins.splt[[i1]]$counts
		    }
		    mu <- mean(qcounts)
		    variance <- var(qcounts)
		    if (names(bins.splt)[i1] == '0-somy') {
		        distr <- 'dgeom'
		        size <- NA
		    }
		    if (variance <= mu) {
		        distr <- 'dbinom'
            size <- dbinom.size(mu, variance)
            prob <- dbinom.prob(mu, variance)
		    } else {
		        distr <- 'dnbinom'
            size <- dnbinom.size(mu, variance)
            prob <- dnbinom.prob(mu, variance)
		    }
		    distributions[[i1]] <- data.frame(type=distr, size=size, prob=prob, mu=mu, variance=variance)
		}
		distributions <- do.call(rbind, distributions)
		rownames(distributions) <- state.labels
		result$distributions <- distributions
  	## Quality info
		qualityInfo <- list(shannon.entropy=qc.entropy(counts), spikiness=qc.spikiness(counts), complexity=attr(result$bins,'qualityInfo')$complexity, bhattacharyya=qc.bhattacharyya(result))
		result$qualityInfo <- qualityInfo
  	## Issue warnings
  	result$warnings <- warlist

  	## Return results
  	return(result)
  
  
}