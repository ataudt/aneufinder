

#' Measures for Karyotype Heterogeneity
#' 
#' Computes measures for karyotype heterogeneity. See the Details section for how these measures are defined.
#'
#' We define \eqn{x} as the vector of copy number states for each position. The number of HMMs is \eqn{S}. The measures are computed for each bin as follows:
#' \describe{
#' \item{Aneuploidy:}{\eqn{D = mean( abs(x-P) )}, where P is the physiological number of chromosomes at that position.}
#' \item{Heterogeneity:}{\eqn{H = sum( table(x) * 0:(length(table(x))-1) ) / S}}
#' }
#'
#' @param hmms A list with \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param normalChromosomeNumbers A named integer vector with physiological copy numbers. This is useful to specify male and female samples, e.g. \code{c('X'=2)} for female samples and \code{c('X'=1,'Y'=1)} for male samples. The assumed default is '2' for all chromosomes.
#' @return A \code{list} with two \code{data.frame}s, containing the karyotype measures $genomewide and $per.chromosome.
#' @author Aaron Taudt
#' @importFrom stats weighted.mean
#' @export
#'@examples
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata/primary-lung/hmms", package="AneuFinder")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata/metastasis-liver/hmms", package="AneuFinder")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
## Compare karyotype measures between the two cancers
#'normal.chrom.numbers <- rep(2, 23)
#'names(normal.chrom.numbers) <- c(1:22,'X')
#'lung <- karyotypeMeasures(lung.files, normalChromosomeNumbers=normal.chrom.numbers)
#'liver <- karyotypeMeasures(liver.files, normalChromosomeNumbers=normal.chrom.numbers)
#'print(lung$genomewide)
#'print(liver$genomewide)
#'
karyotypeMeasures <- function(hmms, normalChromosomeNumbers=NULL) {

	## Load the files
	hmms <- loadHmmsFromFiles(hmms)

	## If all binsizes are the same the consensus template can be chosen equal to the bins
	ptm <- startTimedMessage("Making consensus template ...")
	binsizes <- unlist(lapply(hmms, function(x) { width(x$bins)[1] }))
	# Filter out HMMs where segments of bins$state are NULL
	mask <- !sapply(hmms, function(hmm) { is.null(hmm$segments) | is.null(hmm$bins$state) })
	hmms <- hmms[mask]
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
		consensus <- disjoin(unlist(grlred, use.names=FALSE))
		constates <- matrix(NA, ncol=length(hmms), nrow=length(consensus))
		for (i1 in 1:length(grlred)) {
			grred <- grlred[[i1]]
			splt <- split(grred, mcols(grred)$state)
			multiplicity <- initializeStates(names(splt))$multiplicity
			mind <- as.matrix(findOverlaps(consensus, splt, select='first'))
			constates[,i1] <- multiplicity[names(splt)[mind]]
		}
	}
	stopTimedMessage(ptm)
	

	### Karyotype measures ###
	result <- list()
	S <- ncol(constates)
	## Genomewide
	physioState <- rep(2,length(consensus))
	names(physioState) <- seqnames(consensus)
	if (!is.null(normalChromosomeNumbers)) {
		physioState[names(physioState) %in% names(normalChromosomeNumbers)] <- normalChromosomeNumbers[names(physioState)[names(physioState) %in% names(normalChromosomeNumbers)]]
	}
	consensus$Aneuploidy <- rowMeans(abs(constates - physioState))
	tabs <- apply(constates, 1, function(x) { sort(table(x), decreasing=TRUE) })
	if (is.list(tabs) | is.vector(tabs)) {
		consensus$Heterogeneity <- unlist(lapply(tabs, function(x) { sum(x * 0:(length(x)-1)) })) / S
	} else if (is.matrix(tabs)) {
		consensus$Heterogeneity <- colSums( tabs * 0:(nrow(tabs)-1) ) / S
	}
	weights <- as.numeric(width(consensus))
	result[['genomewide']] <- data.frame(Aneuploidy = stats::weighted.mean(consensus$Aneuploidy, weights),
																			Heterogeneity = stats::weighted.mean(consensus$Heterogeneity, weights))
	## Chromosomes
	consensus.split <- split(consensus, seqnames(consensus))
	weights.split <- split(weights, seqnames(consensus))
	result[['per.chromosome']] <- data.frame(Aneuploidy = unlist(mapply(function(x,y) { stats::weighted.mean(x$Aneuploidy, y) }, consensus.split, weights.split)),
																						Heterogeneity = unlist(mapply(function(x,y) { stats::weighted.mean(x$Heterogeneity, y) }, consensus.split, weights.split)))

	return(result)

}


