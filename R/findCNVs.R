#' Find copy number variations
#'
#' \code{findCNVs} classifies the binned read counts into several states which represent copy-number-variation.
#'
#' \code{findCNVs} uses a 6-state Hidden Markov Model to classify the binned read counts: state 'nullsomy' with a delta function as emission densitiy (only zero read counts), 'monosomy','disomy','trisomy','tetrasomy' and 'multisomy' with negative binomials (see \code{\link{dnbinom}}) as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method. TODO: insert paper
#' @author Aaron Taudt
#' @param method Currently only \code{c('univariate')}.
#' @inheritParams univariate.findCNVs
#' @examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata/BB140820_I_002.bam", package="aneufinder")
#'## Bin the BAM file into bin size 200000bp
#'binned.data <- bam2binned(bamfile, binsize=200000, chromosomes=c(1:22,'X','Y'), GC.correction=FALSE,
#'                          save.as.RData=FALSE)
#'## Fit the Hidden Markov Model
#'model <- findCNVs(binned.data, ID=basename(bamfile), eps=0.1, max.time=60)
#'## Check the fit
#'plot(model, type='histogram')
#' @export
findCNVs <- function(binned.data, ID, method='univariate', eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=10, eps.try=10*eps, num.threads=1, read.cutoff.quantile=0.999, GC.correction=TRUE, strand='*', states=c("zero-inflation","monosomy","disomy","trisomy","tetrasomy","multisomy")) {

	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	message("Find CNVs for ID = ",ID, ":")

	if (method=='univariate') {
		model <- univariate.findCNVs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, GC.correction=GC.correction, strand=strand, states=states)
	}

	attr(model, 'call') <- call
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(model)

}


#' Initialize state factor levels and distributions
#'
#' Initialize the state factor levels and distributions for the specified states.
#'
#' @param states A subset of \code{c("zero-inflation","monosomy","disomy","trisomy","tetrasomy","multisomy")}.
initializeStates <- function(states) {

	possible.states <- c("zero-inflation","monosomy","disomy","trisomy","tetrasomy","multisomy")
	possible.distributions <- factor(c("zero-inflation"='delta',
																			"monosomy"='dnbinom',
																			"disomy"='dnbinom',
																			"trisomy"='dnbinom',
																			"tetrasomy"='dnbinom',
																			"multisomy"='dnbinom'), levels=c('delta','dgeom','dnbinom','dbinom'))
	multiplicity <- c("zero-inflation"=0,
										"monosomy"=1,
										"disomy"=2,
										"trisomy"=3,
										"tetrasomy"=4,
										"multisomy"=5)
	if (any(!(states %in% possible.states))) {
		stop('argument \'states\' accepts only entries from c("zero-inflation","monosomy","disomy","trisomy","tetrasomy","multisomy")')
	}
	labels <- factor(states, levels=possible.states[possible.states %in% states])
	distributions <- possible.distributions[states]
	
	# Return list
	l <- list(labels=labels, distributions=distributions, multiplicity=multiplicity)
	return(l)
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
#'		\item{\code{standard}}{The negative binomial of state 'disomy' will be initialized with \code{mean=mean(reads)}, \code{var=var(reads)}. This procedure usually gives good convergence.}
#'		\item{\code{random}}{Mean and variance of the negative binomial of state 'disomy' will be initialized with random values (in certain boundaries, see source code). Try this if the \code{standard} procedure fails to produce a good fit.}
#'	}
#' @param max.time The maximum running time in seconds for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.
#' @param max.iter The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.
#' @param num.trials The number of trials to find a fit where state 'disomic' is most frequent. Each time, the HMM is seeded with different random initial values.
#' @param eps.try If code num.trials is set to greater than 1, \code{eps.try} is used for the trial runs. If unset, \code{eps} is used.
#' @param num.threads Number of threads to use. Setting this to >1 may give increased performance.
#' @param read.cutoff.quantile A quantile between 0 and 1. Should be near 1. Read counts above this quantile will be set to the read count specified by this quantile. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. Set \code{read.cutoff.quantile=1} in this case.
#' @param GC.correction Either \code{TRUE} or \code{FALSE}. If \code{GC.correction=TRUE}, the GC corrected reads have to be present in the input \code{binned.data}, otherwise a warning is thrown and no GC correction is done.
#' @param strand Run the HMM only for the specified strand. One of \code{c('+', '-', '*')}.
#' @param states A subset or all of \code{c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy")}. This vector defines the states that are used in the Hidden Markov Model. The order of the entries should not be changed.
univariate.findCNVs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, GC.correction=TRUE, strand='*', states=c("zero-inflation","monosomy","disomy","trisomy","tetrasomy","multisomy")) {

	### Define cleanup behaviour ###
	on.exit(.C("R_univariate_cleanup"))

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
	if (check.logical(GC.correction)!=0) stop("argument 'GC.correction' expects a logical (TRUE or FALSE)")
	if (check.strand(strand)!=0) stop("argument 'strand' expects either '+', '-' or '*'")

	warlist <- list()
	if (is.null(eps.try) | num.trials==1) eps.try <- eps

	## Assign variables
	temp <- initializeStates(states)
	state.labels <- temp$labels
	state.distributions <- temp$distributions
	dependent.states.mask <- state.labels %in% c("monosomy","disomy","trisomy","tetrasomy","multisomy")
	numstates <- length(states)
	numbins <- length(binned.data)
	iniproc <- which(init==c("standard","random")) # transform to int
	if (strand=='+') {
		select <- 'preads'
	} else if (strand=='-') {
		select <- 'mreads'
	} else if (strand=='*') {
		select <- 'reads'
	}
	if (GC.correction) {
		if (!(paste0(select,'.gc') %in% names(mcols(binned.data)))) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Cannot use GC-corrected reads because they are not in the binned data. Continuing without GC-correction."))
		} else if (any(is.na(mcols(binned.data)[,paste0(select,'.gc')]))) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Cannot use GC-corrected reads because there are NAs. GC-correction may not be reliable. Continuing without GC-correction."))
		} else if (any(mcols(binned.data)[,paste0(select,'.gc')]<0)) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Cannot use GC-corrected reads because there are negative values. GC-correction is not reliable. Continuing without GC-correction."))
		} else {
			select <- paste0(select,'.gc')
		}
	}
	reads <- mcols(binned.data)[,select]

	### Make return object
		result <- list()
		class(result) <- class.univariate.hmm
		result$ID <- ID
		result$bins <- binned.data
	## Quality info
		qualityInfo <- list(shannon.entropy=qc.entropy(reads), spikyness=qc.spikyness(reads), complexity=attr(result$bins, 'complexity.preseqR'))
		result$qualityInfo <- qualityInfo

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(reads!=0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All reads in data are zero. No HMM done."))
		result$warnings <- warlist
		return(result)
	} else if (any(reads<0)) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some reads in data are negative. No HMM done."))
		result$warnings <- warlist
		return(result)
	}
		

	# Filter high reads out, makes HMM faster
	read.cutoff <- quantile(reads, read.cutoff.quantile)
	names.read.cutoff <- names(read.cutoff)
	read.cutoff <- ceiling(read.cutoff)
	mask <- reads > read.cutoff
	reads[mask] <- read.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		message(paste0("Replaced read counts > ",read.cutoff," (",names.read.cutoff," quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'read.cutoff.quantile=1' to disable this filtering. This filtering was done to increase the speed of the HMM and should not affect the results.\n"))
	}
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		message(paste0("Trial ",i_try," / ",num.trials))

		## Initial parameters
		if (init == 'random') {
			A.initial <- matrix(runif(numstates^2), ncol=numstates)
			A.initial <- sweep(A.initial, 1, rowSums(A.initial), "/")			
			proba.initial <- runif(numstates)
			# Distributions for dependent states
			size.initial <- runif(1, min=0, max=100) * cumsum(dependent.states.mask)
			prob.initial <- runif(1) * dependent.states.mask
		} else if (init == 'standard') {
			A.initial <- matrix(NA, ncol=numstates, nrow=numstates)
			for (irow in 1:numstates) {
				for (icol in 1:numstates) {
					if (irow==icol) { A.initial[irow,icol] <- 0.9 }
					else { A.initial[irow,icol] <- 0.1/(numstates-1) }
				}
			}
			proba.initial <- rep(1/numstates, numstates)
			mean.initial.monosomy <- mean(reads[reads>0])/2
			var.initial.monosomy <- var(reads[reads>0])/2
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
		}
	
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
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
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)

		hmm$eps <- eps.try
		if (num.trials > 1) {
			if (hmm$loglik.delta > hmm$eps) {
				warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": HMM did not converge in trial run ",i_try,"!\n"))
			}
			# Store model in list
			hmm$reads <- NULL
			modellist[[i_try]] <- hmm
			# Check if disomic state is most frequent and stop trials
			if ("disomy" %in% state.labels) {
				idisomy <- which(state.labels=='disomy')
				if (which.max(hmm$weights)==idisomy) {
					break
				}
			}
			init <- 'random'
		} else if (num.trials == 1) {
			if (hmm$loglik.delta > eps) {
				warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": HMM did not converge!\n"))
			}
		}
	}

	if (num.trials > 1) {

		# Select fit with highest weight in state disomic
		# Mathematically we should select the fit with highest loglikelihood. If we think the fit with the highest loglikelihood is incorrect, we should change the underlying model. However, this is very complex and we choose to select a fit that we think is (more) correct, although it has not the highest support given our (imperfect) model.
		if ("disomy" %in% state.labels) {
			idisomy <- which(state.labels=='disomy')
			indexmax <- which.max(unlist(lapply(lapply(modellist,'[[','weights'), '[', idisomy)))
		} else {
			indexmax <- which.max(unlist(lapply(modellist,'[[','loglik'))) # fit with highest loglikelihood
		}
		hmm <- modellist[[indexmax]]

		# Check if size and prob parameter are correct
		if (any(is.na(hmm$size) | is.nan(hmm$size) | is.infinite(hmm$size) | is.na(hmm$prob) | is.nan(hmm$prob) | is.infinite(hmm$prob))) {
			hmm$error <- 3
		} else {

			# Rerun the HMM with different epsilon and initial parameters from trial run
			message(paste0("Rerunning trial ",indexmax," with eps = ",eps))
			hmm <- .C("R_univariate_hmm",
				reads = as.integer(reads), # int* O
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
				read.cutoff = as.integer(read.cutoff) # int* read_cutoff
			)
		}

	}

	### Make return object ###
	## Check for errors
		if (hmm$error == 0) {
		## Bin coordinates and states ###
			result$bins$state <- state.labels[hmm$states]
		## Segmentation
			message("Making segmentation ...", appendLF=F)
			ptm <- proc.time()
			gr <- result$bins
			red.gr.list <- GRangesList()
			for (state in state.labels) {
				red.gr <- GenomicRanges::reduce(gr[gr$state==state])
				mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels)),length(red.gr))
				if (length(red.gr) > 0) {
					red.gr.list[[length(red.gr.list)+1]] <- red.gr
				}
			}
			red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
			result$segments <- red.gr
			seqlengths(result$segments) <- seqlengths(binned.data)
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
			names(result$startProbs) <- paste0("P(",state.labels,")")
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- paste0("P(",state.labels,")")
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
				result$distributions.initial <- distributions
		## Convergence info
			convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec, error=hmm$error)
			result$convergenceInfo <- convergenceInfo
		## GC correction
			result$GC.correction <- grepl('gc', select)
		## Quality info
			qualityInfo <- list(shannon.entropy=qc.entropy(reads), spikyness=qc.spikyness(reads), complexity=attr(result$bins, 'complexity.preseqR'))
			result$qualityInfo <- qualityInfo
		} else if (hmm$error == 1) {
			warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'read.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
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


