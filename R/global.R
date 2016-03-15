

#' @useDynLib AneuFinder, .registration = TRUE, .fixes = ""
#' @import IRanges
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import AneuFinderData
NULL

# =======================================================
# Some global variables that can be used in all functions
# =======================================================
## Class names
class.univariate.hmm <- "aneuHMM"
class.multivariate.hmm <- "aneuMultiHMM"
class.bivariate.hmm <- "aneuBiHMM"
class.hmm.list <- "aneuHMM.list"

