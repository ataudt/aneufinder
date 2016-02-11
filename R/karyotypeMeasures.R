

#' Measures for Karyotype Heterogeneity
#' 
#' Computes measures for karyotype heterogeneity. See the Details section for how these measures are defined.
#'
#' We define \eqn{x} as the vector of copy number states for each position and HMM. The number of HMMs is \eqn{S}.
#' \describe{
#' \item{Divergence from disomic:}{\eqn{D = }}
#' \item{Simple karyotype heterogeneity:}{\eqn{H = table(x) * 0:(length(table(x))-1)}}
#' \item{Entropic karyotype heterogeneity:}{\eqn{H = S*log(S) - S + sum(table(x)-table(x)*log(table(x)))}}
#' }
#'
#' @param hmms A list with \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param normalChromosomeNumbers A named integer vector with physiological copy numbers. This is useful to specify male and female samples, e.g. \code{c('X'=2)} for female samples and \code{c('X'=1,'Y'=1)} for male samples. The assumed default is '2' for all chromosomes.
#' @author Aaron Taudt
#' @export
karyotypeMeasures <- function(hmms, normalChromosomeNumbers=NULL) {

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)

	## If all binsizes are the same the consensus template can be chosen equal to the bins
	message("Making consensus template ...", appendLF=F); ptm <- proc.time()
	binsizes <- unlist(lapply(hmms, function(x) { width(x$bins)[1] }))
	if (all(binsizes==binsizes[1])) {
		consensus <- hmms[[1]]$bins
		mcols(consensus) <- NULL
		constates <- matrix(NA, ncol=length(hmms), nrow=length(consensus))
		for (i1 in 1:length(hmms)) {
			hmm <- hmms[[i1]]
			multiplicity <- initializeStates(levels(hmm$bins$state))$multiplicity
			constates[,i1] <- multiplicity[as.character(hmm$bins$state)]
		}
	} else { # binsizes differ
		## Get segments from list
		grlred <- GRangesList()
		for (hmm in hmms) {
			if (!is.null(hmm$segments)) {
				grlred[[hmm$ID]] <- hmm$segments
			}
		}
		## Consensus template
		consensus <- disjoin(unlist(grlred))
		constates <- matrix(NA, ncol=length(hmms), nrow=length(consensus))
		for (i1 in 1:length(grlred)) {
			grred <- grlred[[i1]]
			splt <- split(grred, mcols(grred)$state)
			multiplicity <- initializeStates(names(splt))$multiplicity
			mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
			constates[,i1] <- multiplicity[names(splt)[mind]]
		}
	}
	meanstates <- apply(constates, 1, mean, na.rm=T)
	mcols(consensus)$meanstate <- meanstates
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	

	### Karyotype measures ###
	result <- list()
	S <- ncol(constates)
	## Genomewide
	physioState <- rep(2,length(consensus))
	names(physioState) <- seqnames(consensus)
	if (!is.null(normalChromosomeNumbers)) {
		physioState[names(physioState) %in% names(normalChromosomeNumbers)] <- normalChromosomeNumbers[names(physioState)[names(physioState) %in% names(normalChromosomeNumbers)]]
	}
	consensus$divergenceFromDisomic <- abs(consensus$meanstate - physioState)
	tabs <- apply(constates, 1, function(x) { sort(table(x), decreasing=T) })
	consensus$simpleHeterogeneity <- unlist(lapply(tabs, function(x) { sum(x * 0:(length(x)-1)) })) / S
	consensus$entropicHeterogeneity <- S*log(S) - S + unlist(lapply(tabs, function(x) { sum(x-x*log(x)) }))
	weights <- as.numeric(width(consensus))
	result[['genomewide']] <- data.frame(divergenceFromDisomic = weighted.mean(consensus$divergenceFromDisomic, weights),
																			simpleHeterogeneity = weighted.mean(consensus$simpleHeterogeneity, weights),
																			entropicHeterogeneity = weighted.mean(consensus$entropicHeterogeneity, weights))
	## Chromosomes
	consensus.split <- split(consensus, seqnames(consensus))
	weights.split <- split(weights, seqnames(consensus))
	result[['per.chromosome']] <- data.frame(divergenceFromDisomic = unlist(mapply(function(x,y) { weighted.mean(x$divergenceFromDisomic, y) }, consensus.split, weights.split)),
																						simpleHeterogeneity = unlist(mapply(function(x,y) { weighted.mean(x$simpleHeterogeneity, y) }, consensus.split, weights.split)),
																						entropicHeterogeneity = unlist(mapply(function(x,y) { weighted.mean(x$entropicHeterogeneity, y) }, consensus.split, weights.split)))

	return(result)

}


