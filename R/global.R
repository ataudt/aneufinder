#' @useDynLib aneufinder
#' @import GenomicRanges
#' @import IRanges
NULL

# =======================================================
# Some global variables that can be used in all functions
# =======================================================
class.univariate.hmm <- "aneuHMM"
class.multivariate.hmm <- "aneuMultiHMM"
class.bivariate.hmm <- "aneuBiHMM"
class.hmm.list <- "aneufinder.hmm.list"

state.labels <- factor(c("nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"), levels=c("nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"))
dependent.states.mask <- state.labels %in% c("monosomy","disomy","trisomy","tetrasomy","multisomy")
state.distributions <- factor(c('delta','dnbinom','dnbinom','dnbinom','dnbinom','dnbinom'), levels=c('delta','dgeom','dnbinom','dbinom'))
state.colors <- c("mapped"="gray68","nullsomy"="gray90","monosomy"="gold3","disomy"="springgreen3","trisomy"="orangered1","tetrasomy"="orangered4","multisomy"="purple3","total"="black")
state.colors <- c(mapped='#ADADAD',
									"zero-inflation"='#E5E5E5',
									nullsomy='#E5E5E5',
									monosomy='#A3AD00',
									disomy='#078735',
									trisomy='#DE3000',
									tetrasomy='#660602',
									multisomy='#0B62AB',
									total='#000000')
get.state.labels <- function() { return(state.labels) }
get.state.colors <- function() { return(state.colors[as.character(state.labels)]) }
 
state.labels.SCE <- factor(c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"), levels=c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy"))
multiplicity.SCE <- c(0,0,1,2,3,4,5)
names(multiplicity.SCE) <- state.labels.SCE
dependent.states.mask.SCE <- state.labels.SCE %in% c("monosomy","disomy","trisomy","tetrasomy","multisomy")
state.distributions.SCE <- factor(c('delta','dgeom','dnbinom','dnbinom','dnbinom','dnbinom','dnbinom'), levels=c('delta','dgeom','dnbinom','dbinom'))
get.state.labels.SCE <- function() { return(state.labels.SCE) }

# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
dnbinom.size <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

dnbinom.prob <- function(mean, variance) {
	return(mean/variance)
}

dnbinom.mean <- function(size, prob) {
	return(size/prob - size)
}

dnbinom.variance <- function(size, prob) {
	return( (size - prob*size) / prob^2 )
}

dgeom.mean <- function(prob) {
	return( (1-prob)/prob )
}

dgeom.variance <- function(prob) {
	return( (1-prob)/prob^2 )
}

dbinom.size <- function(mean, variance) {
	return( mean^2/(mean-variance) )
}

dbinom.prob <- function(mean, variance) {
	return( (mean-variance)/mean )
}

dbinom.mean <- function(size, prob) {
	return( size*prob )
}

dbinom.variance <- function(size, prob) {
	return( size*prob * (1-prob) )
}


