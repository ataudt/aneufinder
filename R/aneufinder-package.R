

#' Copy-number detection in WGSCS and Strand-Seq data
#'
#' CNV detection in whole-genome single cell sequencing (WGSCS) and Strand-seq data using a Hidden Markov Model. The package implements CNV detection, commonly used plotting functions, export to BED format for upload to genome browsers, and measures for assessment of karyotype heterogeneity and quality metrics.
#'
#' The main function of this package is \code{\link{Aneufinder}} and produces several plots and browser files. If you want to have more fine-grained control over the different steps (binning, GC-correction, HMM, plotting) check the vignette \href{../doc/AneuFinder.pdf}{Introduction to AneuFinder}.
#'
#' @author Aaron Taudt, David Porubsky
#' @docType package
#' @name AneuFinder-package
#' @aliases AneuFinder
NULL
