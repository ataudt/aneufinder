

#' Copy-number detection in whole-genome single-cell sequencing data
#'
#' CNV detection in whole-genome single cell sequencing data using a Hidden Markov Model. The package implements CNV detection, commonly used plotting functions, export to BED format for upload to genome browsers, and measures for assessment of karyotype heterogeneity and quality metrics.
#'
#' The main function of this package is \code{\link{Aneufinder}} and produces several plots and browser files. If you want to have more fine-grained control over the different steps (binning, GC-correction, HMM, plotting) check the vignette \href{../doc/aneufinder.pdf}{Introduction to aneufinder}.
#'
#' @author Aaron Taudt, David Porubsky
#' @docType package
#' @name aneufinder-package
#' @aliases aneufinder
NULL
