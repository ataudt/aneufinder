#' Find sister chromatid exchanges
#'
#' \code{findSCEs} classifies the binned read counts into several states which represent the number of chromatids on each strand.
#'
#' \code{findSCEs} uses a 6-state Hidden Markov Model to classify the binned read counts: state 'nullsomy' with a delta function as emission densitiy (only zero read counts), 'monosomy','disomy','trisomy','tetrasomy' and 'multisomy' with negative binomials (see \code{\link{dnbinom}}) as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method. TODO: insert paper
#' @author Aaron Taudt
#' @inheritParams univariate.findSCEs
#' @return An \code{\link{aneuBiHMM}} object.
#' @examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata/BB140820_I_002.bam", package="aneufinder")
#'## Bin the BAM file into bin size 200000bp
#'binned.data <- bam2binned(bamfile, binsize=200000, chromosomes=c(1:22,'X','Y'), GC.correction=FALSE,
#'                          save.as.RData=FALSE)
#'## Fit the Hidden Markov Model
#'model <- findSCEs(binned.data, ID=basename(bamfile), eps=0.1, max.time=60)
#'## Check the fit
#'plot(model, type='histogram')
#' @export
findSCEs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=10, eps.try=10*eps, num.threads=1, read.cutoff.quantile=0.999, GC.correction=TRUE, strand='*', allow.odd.states=TRUE) {

	call <- match.call()
	underline <- paste0(rep('=',sum(nchar(call[[1]]))+3), collapse='')
	message("\n",call[[1]],"():")
	message(underline)
	ptm <- proc.time()
	message("Find CNVs for ID = ",ID, ":")

	model <- bivariate.findSCEs(binned.data, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=read.cutoff.quantile, GC.correction=GC.correction, allow.odd.states=allow.odd.states)

	attr(model, 'call') <- call
	time <- proc.time() - ptm
	message("Time spent in ", call[[1]],"(): ",round(time[3],2),"s")
	return(model)

}


#' Find copy number variations (univariate)
#'
#' \code{findSCEs} classifies the input binned read counts into several states which represent copy-number-variation.
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
#' @return An \code{\link{aneuHMM}} object.
univariate.findSCEs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, GC.correction=TRUE, strand='*') {

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
	if (is.null(eps.try)) eps.try <- eps

	## Assign variables
# 	state.labels.SCE # assigned in global.R
	numstates <- length(state.labels.SCE)
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
			size.initial <- runif(1, min=0, max=100) * cumsum(dependent.states.mask.SCE)
			prob.initial <- runif(1) * dependent.states.mask.SCE
		} else if (init == 'standard') {
			A.initial <- matrix(NA, ncol=numstates, nrow=numstates)
			for (irow in 1:numstates) {
				for (icol in 1:numstates) {
					if (irow==icol) { A.initial[irow,icol] <- 0.9 }
					else { A.initial[irow,icol] <- 0.1/(numstates-1) }
				}
			}
			proba.initial <- rep(1/numstates, numstates)
			mean.initial <- mean(reads[reads>0]) * cumsum(dependent.states.mask.SCE)
			var.initial <- var(reads[reads>0]) * cumsum(dependent.states.mask.SCE)
			size.initial <- rep(0,numstates)
			prob.initial <- rep(0,numstates)
			mask <- dependent.states.mask.SCE
			size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
			prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
		}
		# Assign initials for the nullsomy distribution
		size.initial[2] <- 1
		prob.initial[2] <- 0.5
	
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			state.labels.SCE = as.integer(state.labels.SCE), # int* state.labels.SCE
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
			distr.type = as.integer(state.distributions.SCE), # int* distr_type
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
			# Check if monosomic is more frequent than trisomic
			imono <- which(state.labels.SCE=='monosomy')
			itri <- which(state.labels.SCE=='trisomy')
			if (hmm$weights[imono]>hmm$weights[itri]) {
				break
			}
			init <- 'random'
		} else if (num.trials == 1) {
			if (hmm$loglik.delta > eps) {
				warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": HMM did not converge!\n"))
			}
		}
	}

	if (num.trials > 1) {

		# Select fit with highest weight in state disomy
		# Mathematically we should select the fit with highest loglikelihood. If we think the fit with the highest loglikelihood is incorrect, we should change the underlying model. However, this is very complex and we choose to select a fit that we think is (more) correct, although it has not the highest support given our (imperfect) model.
		istate <- which(state.labels.SCE=='disomy')
		indexmax <- which.max(unlist(lapply(lapply(modellist,'[[','weights'), '[', istate)))
		hmm <- modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		message(paste0("Rerunning trial ",indexmax," with eps = ",eps))
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # int* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			state.labels.SCE = as.integer(state.labels.SCE), # int* state.labels.SCE
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
			distr.type = as.integer(state.distributions.SCE), # int* distr_type
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

	### Make return object ###
	## Check for errors
		if (hmm$error == 0) {
		## Bin coordinates and states ###
			result$bins$state <- state.labels.SCE[hmm$states]
		## Segmentation
			message("Making segmentation ...", appendLF=F)
			ptm <- proc.time()
			gr <- result$bins
			red.gr.list <- GRangesList()
			for (state in state.labels.SCE) {
				red.gr <- GenomicRanges::reduce(gr[gr$state==state])
				if (length(red.gr) > 0) {
					mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels.SCE)),length(red.gr))
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
			names(result$weights) <- state.labels.SCE
			# Transition matrices
			transitionProbs <- matrix(hmm$A, ncol=hmm$num.states)
			rownames(transitionProbs) <- state.labels.SCE
			colnames(transitionProbs) <- state.labels.SCE
			result$transitionProbs <- transitionProbs
			transitionProbs.initial <- matrix(hmm$A.initial, ncol=hmm$num.states)
			rownames(transitionProbs.initial) <- state.labels.SCE
			colnames(transitionProbs.initial) <- state.labels.SCE
			result$transitionProbs.initial <- transitionProbs.initial
			# Initial probs
			result$startProbs <- hmm$proba
			names(result$startProbs) <- paste0("P(",state.labels.SCE,")")
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- paste0("P(",state.labels.SCE,")")
			# Distributions
				distributions <- data.frame()
				distributions.initial <- data.frame()
				for (idistr in 1:length(hmm$distr.type)) {
					distr <- levels(state.distributions.SCE)[hmm$distr.type[idistr]]
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
				rownames(distributions) <- state.labels.SCE
				rownames(distributions.initial) <- state.labels.SCE
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
		}

	## Issue warnings
	result$warnings <- warlist

	## Return results
	return(result)
}


bivariate.findSCEs <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, GC.correction=TRUE, allow.odd.states=TRUE) {

	## Debugging
# 	ID = 'test'
# 	eps = 1
# 	init='standard'
# 	max.time=-1
# 	max.iter=-1
# 	num.trials=1
# 	eps.try=NULL
# 	num.threads=1
# 	read.cutoff.quantile=0.999
# 	GC.correction=FALSE
#   allow.odd.states=T

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

	warlist <- list()
	if (is.null(eps.try)) eps.try <- eps

	## Variables
	num.bins <- length(binned.data)

	## Get reads
	select <- 'reads'
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
	reads <- matrix(c(mcols(binned.data)[,paste0('m',select)], mcols(binned.data)[,paste0('p',select)]), ncol=2, dimnames=list(bin=1:num.bins, strand=c('minus','plus')))
	maxreads <- max(reads)

	## Filter high reads out, makes HMM faster
	read.cutoff <- quantile(reads, read.cutoff.quantile)
	names.read.cutoff <- names(read.cutoff)
	read.cutoff <- ceiling(read.cutoff)
	mask <- reads > read.cutoff
	reads[mask] <- read.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		message(paste0("Replaced read counts > ",read.cutoff," (",names.read.cutoff," quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'read.cutoff.quantile=1' to disable this filtering. This filtering was done to increase the speed of the HMM and should not affect the results.\n"))
	}
	
	### Make return object
		result <- list()
		class(result) <- class.bivariate.hmm
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

	### Stack the strands and run one univariate findSCEs
	message("")
	message(paste(rep('-',getOption('width')), collapse=''))
	message("Running univariate")
	binned.data.minus <- binned.data
	strand(binned.data.minus) <- '-'
	binned.data.minus$reads <- binned.data.minus$mreads
	binned.data.minus$reads.gc <- binned.data.minus$mreads.gc
	binned.data.plus <- binned.data
	strand(binned.data.plus) <- '+'
	binned.data.plus$reads <- binned.data.plus$preads
	binned.data.plus$reads.gc <- binned.data.plus$preads.gc
	binned.data.stacked <- c(binned.data.minus, binned.data.plus)
	mask.attributes <- c(grep('complexity', names(attributes(binned.data)), value=T), 'spikyness', 'shannon.entropy')
	attributes(binned.data.stacked)[mask.attributes] <- attributes(binned.data)[mask.attributes]

	model.stacked <- univariate.findSCEs(binned.data.stacked, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff.quantile=1, GC.correction=GC.correction)
	model.minus <- model.stacked
	model.minus$bins <- model.minus$bins[strand(model.minus$bins)=='-']
	model.minus$segments <- model.minus$segments[strand(model.minus$segments)=='-']
	model.plus <- model.stacked
	model.plus$bins <- model.plus$bins[strand(model.plus$bins)=='+']
	model.plus$segments <- model.plus$segments[strand(model.plus$segments)=='+']

	models <- list(minus=model.minus, plus=model.plus)

	### Prepare the multivariate HMM
	message("")
	message(paste(rep('-',getOption('width')), collapse=''))
	message("Preparing bivariate HMM\n")

		## Extract reads and other stuff
		distributions <- lapply(models, '[[', 'distributions')
		weights <- lapply(models, '[[', 'weights')
		num.uni.states <- length(models[[1]]$weights)
		num.models <- length(models)
		num.bins <- length(models[[1]]$bins)
		comb.states <- vector()
		levels.state <- levels(models[[1]]$bins$state)
		for (i1 in 1:length(levels.state)) {
			for (i2 in 1:length(levels.state)) {
				comb.state <- paste(levels.state[i1], levels.state[i2])
				if (allow.odd.states) {
					comb.states[length(comb.states)+1] <- comb.state
				} else {
					state.multiplicity <- multiplicity.SCE[levels.state[i1]] + multiplicity.SCE[levels.state[i2]]
					# Only even states and null-mono states (for sex chromosomes)
					if (state.multiplicity %% 2 == 0 | state.multiplicity == 1) {
						comb.states[length(comb.states)+1] <- comb.state
					}
				}
			}
		}
		comb.states <- factor(comb.states, levels=comb.states)
		states.list <- list()
		for (model in models) { states.list[[length(states.list)+1]] <- model$bins$state }
		comb.states.per.bin <- factor(do.call(paste, states.list), levels=levels(comb.states))
		num.comb.states <- length(comb.states)

		## Pre-compute z-values for each number of reads
		message("Computing pre z-matrix...", appendLF=F)
		z.per.read <- array(NA, dim=c(maxreads+1, num.models, num.uni.states), dimnames=list(reads=0:maxreads, strand=names(models), state=levels(models[[1]]$bins$state)))
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
				} else if (model$distributions[istate,'type']=='dgeom') {
					prob <- model$distributions[istate,'prob']
					u <- pgeom(xreads, prob)
				}
				qnorm_u <- qnorm(u)
				mask <- qnorm_u==Inf
				qnorm_u[mask] <- qnorm(1-1e-16)
				z.per.read[ , istrand, istate] <- qnorm_u
			}
		}
		message(" done")

		## Compute the z matrix
		message("Transfering values into z-matrix...", appendLF=F)
		z.per.bin = array(NA, dim=c(num.bins, num.models, num.uni.states), dimnames=list(bin=1:num.bins, strand=names(models), state=levels(models[[1]]$bins$state)))
		for (istrand in 1:num.models) {
			model <- models[[istrand]]
			for (istate in 1:num.uni.states) {
				z.per.bin[ , istrand, istate] = z.per.read[reads[,istrand]+1, istrand, istate]
			}
		}
		remove(z.per.read)
		message(" done")

		## Calculate correlation matrix
		message("Computing inverse of correlation matrix...", appendLF=F)
		correlationMatrix <- array(0, dim=c(num.models,num.models,num.comb.states), dimnames=list(strand=names(models), strand=names(models), comb.state=comb.states))
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
			}, warning = function(war) {
				correlationMatrix[,,comb.state] <<- diag(num.models)
				determinant[comb.state] <<- det( correlationMatrix[,,comb.state] )
				correlationMatrixInverse[,,comb.state] <<- solve(correlationMatrix[,,comb.state])
				usestateTF[comb.state] <<- TRUE
				war
			}, error = function(err) {
				correlationMatrix[,,comb.state] <<- diag(num.models)
				determinant[comb.state] <<- det( correlationMatrix[,,comb.state] )
				correlationMatrixInverse[,,comb.state] <<- solve(correlationMatrix[,,comb.state])
				usestateTF[comb.state] <<- TRUE
				err
			})
		}
		message(" done")

		# Use only states with valid correlationMatrixInverse
		correlationMatrix = correlationMatrix[,,usestateTF]
		correlationMatrixInverse = correlationMatrixInverse[,,usestateTF]
		comb.states = comb.states[usestateTF]
		comb.states <- droplevels(comb.states)
		determinant = determinant[usestateTF]
		num.comb.states <- length(comb.states)

		## Calculate multivariate densities for each state
		message("Calculating multivariate densities...", appendLF=F)
		densities <- matrix(1, ncol=num.comb.states, nrow=num.bins, dimnames=list(bin=1:num.bins, comb.state=comb.states))
		for (comb.state in comb.states) {
			istate <- which(comb.state==comb.states)
			state <- strsplit(comb.state, ' ')[[1]]
			z.temp <- matrix(NA, ncol=num.models, nrow=num.bins)
			product <- 1
			for (istrand in 1:num.models) {
				z.temp[,istrand] <- z.per.bin[, istrand, state[istrand]]
				if (models[[istrand]]$distributions[state[istrand],'type'] == 'dnbinom') {
					size <- models[[istrand]]$distributions[state[istrand],'size']
					prob <- models[[istrand]]$distributions[state[istrand],'prob']
					product <- product * dnbinom(reads[,istrand], size, prob)
				} else if (models[[istrand]]$distributions[state[istrand],'type'] == 'dgeom') {
					prob <- models[[istrand]]$distributions[state[istrand],'prob']
					product <- product * dgeom(reads[,istrand], prob)
				} else if (models[[istrand]]$distributions[state[istrand],'type'] == 'delta') {
					product <- product * ifelse(reads[,istrand]==0, 1, 0)
				}
			}
			exponent <- -0.5 * apply( ( z.temp %*% (correlationMatrixInverse[ , , istate] - diag(num.models)) ) * z.temp, 1, sum)
			densities[,istate] <- product * determinant[istate]^(-0.5) * exp( exponent )
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
		message(" done\n", appendLF=F)
		
	### Define cleanup behaviour ###
	on.exit(.C("R_multivariate_cleanup", as.integer(num.comb.states)))

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
		war <- warning(paste0("ID = ",ID,": HMM did not converge!\n"))
	}
	if (hmm$error == 1) {
		warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'read.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
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
		matrix.states <- matrix(unlist(strsplit(as.character(result$bins$state), split=' ')), byrow=T, ncol=num.models, dimnames=list(bin=1:num.bins, strand=names(models)))
		names <- c('mstate','pstate')
		for (i1 in 1:num.models) {
			mcols(result$bins)[names[i1]] <- factor(matrix.states[,i1], levels=levels(models[[i1]]$bins$state))
		}
	## Segmentation
		message("Making segmentation ...", appendLF=F)
		ptm <- proc.time()
		gr <- result$bins
		red.gr.list <- GRangesList()
		for (state in comb.states) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			if (length(red.gr)>0) {
				mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
				red.gr.list[[length(red.gr.list)+1]] <- red.gr
			}
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		result$segments <- red.gr
		seqlengths(result$segments) <- seqlengths(result$bins)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
	## CNV state for both strands combined
		# Bins
		state <- multiplicity.SCE[result$bins$mstate] + multiplicity.SCE[result$bins$pstate]
		state[state>max(multiplicity.SCE)] <- max(multiplicity.SCE)
		multiplicity.SCE.inverse <- names(multiplicity.SCE)
		names(multiplicity.SCE.inverse) <- multiplicity.SCE
		state <- multiplicity.SCE.inverse[as.character(state)]
		state[(result$bins$mstate=='nullsomy' | result$bins$pstate=='nullsomy') & state=='zero-inflation'] <- 'nullsomy'
    result$bins$state <- factor(state, levels=state.labels.SCE)
		# Segments
		str <- strsplit(as.character(result$segments$state),' ')
		result$segments$mstate <- factor(unlist(lapply(str, '[[', 1)), levels=state.labels.SCE)
		result$segments$pstate <- factor(unlist(lapply(str, '[[', 2)), levels=state.labels.SCE)
		state <- multiplicity.SCE[result$segments$mstate] + multiplicity.SCE[result$segments$pstate]
		state[state>max(multiplicity.SCE)] <- max(multiplicity.SCE)
		multiplicity.SCE.inverse <- names(multiplicity.SCE)
		names(multiplicity.SCE.inverse) <- multiplicity.SCE
		state <- multiplicity.SCE.inverse[as.character(state)]
		state[(result$segments$mstate=='nullsomy' | result$segments$pstate=='nullsomy') & state=='zero-inflation'] <- 'nullsomy'
    result$segments$state <- factor(state, levels=state.labels.SCE)
		## Parameters
			# Weights
			tstates <- table(result$bins$state)
			result$weights <- tstates/sum(tstates)
			result$weights.univariate <- weights
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
		## Quality info
			qualityInfo <- list(shannon.entropy=qc.entropy(reads), spikyness=qc.spikyness(reads), complexity=attr(result$bins, 'complexity.preseqR'))
			result$qualityInfo <- qualityInfo
	} else if (hmm$error == 1) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library! The following factors are known to cause this error: 1) Your read counts contain very high numbers. Try again with a lower value for 'read.cutoff.quantile'. 2) Your library contains too few reads in each bin. 3) Your library contains reads for a different genome than it was aligned to."))
	} else if (hmm$error == 2) {
		warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": An error occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your library"))
	}

	## Issue warnings
	result$warnings <- warlist

	## Return results
	return(result)
}


#' Filter segments by minimal size
#'
#' \code{filterSegments} filters out segments below a specified minimal segment size. This can be useful to get rid of boundary effects from the Hidden Markov approach.
#'
#' @param model A \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} object.
#' @param min.seg.width The minimum segment width in base-pairs. The default \code{NULL} will be twice the bin size of the given \code{model}.
#' @return The input \code{model} with adjusted segments.
#' @author Aaron Taudt
#' @export
filterSegments <- function(model, min.seg.width=NULL) {
	
	if (is.null(min.seg.width)) {
		min.seg.width <- 2*width(model$bins)[1]
	}

# 	## Reassign state column for segmentation
# 	model$bins$state <- factor(paste(model$bins$mstate, model$bins$pstate))
# 	model$segments$state <- factor(paste(model$segments$mstate, model$segments$pstate))
# 
# 	grl <- GRangesList()
# 	for (chrom in seqlevels(model$segments)) {
# 		gr <- model$segments[seqnames(model$segments)==chrom]
# 		## At the beginning of each chromosome, assign the state after to segments smaller than min.seg.width
# 		index <- which(width(gr) < min.seg.width)
# 		if (length(index)==0) {
# 			grl[[chrom]] <- gr
# 			next
# 		}
# 		post.index <- index+1
# 		# Filter consecutive indices
# 		if (length(index)>1) {
# 			for (i1 in (length(index)-1):1) {
# 				if (index[i1]==index[i1+1]-1) {
# 					post.index[i1] <- post.index[i1+1]
# 				}
# 			}
# 		}
# 		# Only beginning of chromosome
# 		mask <- post.index==post.index[1]
# 		post.index <- post.index[mask]
# 		index <- index[mask]
# 		if (index[1]==1) {
# 			gr$state[index] <- gr$state[post.index]
# 			## Redo segmentation
# 			red.gr.list <- GRangesList()
# 			for (state in levels(gr$state)) {
# 				red.gr <- GenomicRanges::reduce(gr[gr$state==state])
# 				if (length(red.gr)>0) {
# 					mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
# 					red.gr.list[[length(red.gr.list)+1]] <- red.gr
# 				}
# 			}
# 			red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
# 			gr <- red.gr
# 			seqlengths(gr) <- seqlengths(model$bins)
# 		}
# 		## For the rest of each chromosome, assign the state before to segments smaller than min.seg.width
# 		index <- which(width(gr) < min.seg.width)
# 		prev.index <- index-1
# 		# Filter consecutive indices
# 		if (length(index)>1) {
# 			for (i1 in 2:length(index)) {
# 				if (index[i1]==index[i1-1]+1) {
# 					prev.index[i1] <- prev.index[i1-1]
# 				}
# 			}
# 		}
# 		gr$state[index] <- gr$state[prev.index]
# 
# 		grl[[chrom]] <- gr
# 	}
# 		
# 	## Redo segmentation
# 	gr <- unlist(grl)
# 	red.gr.list <- GRangesList()
# 	for (state in levels(gr$state)) {
# 		red.gr <- GenomicRanges::reduce(gr[gr$state==state])
# 		if (length(red.gr)>0) {
# 			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
# 			red.gr.list[[length(red.gr.list)+1]] <- red.gr
# 		}
# 	}
# 	red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
# 	model$segments <- red.gr
# 	seqlengths(model$segments) <- seqlengths(model$bins)
# 
# 	## CNV state for both strands combined
# 	# Bins
# 	state <- multiplicity.SCE[model$bins$mstate] + multiplicity.SCE[model$bins$pstate]
# 	state[state>max(multiplicity.SCE)] <- max(multiplicity.SCE)
# 	multiplicity.SCE.inverse <- names(multiplicity.SCE)
# 	names(multiplicity.SCE.inverse) <- multiplicity.SCE
# 	state <- multiplicity.SCE.inverse[as.character(state)]
# 	state[(model$bins$mstate=='nullsomy' | model$bins$pstate=='nullsomy') & state=='zero-inflation'] <- 'nullsomy'
# 	model$bins$state <- factor(state, levels=state.labels.SCE)
# 	# Segments
# 	str <- strsplit(as.character(model$segments$state),' ')
# 	model$segments$mstate <- factor(unlist(lapply(str, '[[', 1)), levels=state.labels.SCE)
# 	model$segments$pstate <- factor(unlist(lapply(str, '[[', 2)), levels=state.labels.SCE)
# 	state <- multiplicity.SCE[model$segments$mstate] + multiplicity.SCE[model$segments$pstate]
# 	state[state>max(multiplicity.SCE)] <- max(multiplicity.SCE)
# 	multiplicity.SCE.inverse <- names(multiplicity.SCE)
# 	names(multiplicity.SCE.inverse) <- multiplicity.SCE
# 	state <- multiplicity.SCE.inverse[as.character(state)]
# 	state[(model$segments$mstate=='nullsomy' | model$segments$pstate=='nullsomy') & state=='zero-inflation'] <- 'nullsomy'
# 	model$segments$state <- factor(state, levels=state.labels.SCE)

	model$bins <- model$bins[width(model$bins) >= min.seg.width]
	model$segments <- model$segments[width(model$segments) >= min.seg.width]

	return(model)
}

#' Get SCE coordinates
#'
#' Extracts the coordinates of a sister chromatid exchanges (SCE) from an \code{\link{aneuBiHMM}} object.
#'
#' @param model An \code{\link{aneuBiHMM}} object.
#' @return A \code{\link{GRanges}} object containing the SCE coordinates.
#' @author Aaron Taudt
#' @export
getSCEcoordinates <- function(model, resolution=c(1e6)) {

	sce <- GRanges()
	for (ires in resolution) {
		model.filtered <- filterSegments(model, min.seg.width=ires)
		segments <- model.filtered$segments
		if (is.null(segments)) {
			return(sce)
		}
		segments <- segments[segments$state != 'zero-inflation' & segments$state != 'nullsomy']
		segments.split <- split(segments, seqnames(segments))
		for (chrom in names(segments.split)) {
			segments <- segments.split[[chrom]]
			if (length(segments)>1) {
				multiplicity.minus <- multiplicity.SCE[segments$mstate]
				multiplicity.plus <- multiplicity.SCE[segments$pstate]
				multiplicity.both <- multiplicity.SCE[segments$state]
				mask.minus <- c(0,diff(as.vector(multiplicity.minus))) != 0
				mask.both <- c(-1,diff(as.vector(multiplicity.both))) == 0
				mask <- mask.minus & mask.both
	# segments$mask.minus <- mask.minus
	# segments$mask.both <- mask.both
	# segments$mask <- mask
	# segments$w <- width(segments)
				if (length(which(mask))>0) {
					sce <- c(sce, segments[mask])
				}
			}
		}
	}
	mcols(sce) <- NULL
	end(sce) <- start(sce)
	sce <- reduce(sce)

	return(sce)
}

