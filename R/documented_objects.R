

#' Binned read counts
#'
#' A \link{GRanges-class} object which contains binned read counts as meta data column \code{reads}. It is output of the various \link{binning} functions.
#' @name binned.data
NULL

#' Hidden Markov Model
#'
#' The \code{aneuHMM} object is output of the function \code{\link{findCNVs}} and is basically a list with various entries. The class() attribute of this list was set to "aneuHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
#'
#' @name aneuHMM
#' @return
#' \item{ID}{An identifier that is used in various \pkg{\link{AneuFinder}} functions.}
#' \item{bins}{
#' A \link{GRanges-class} object containing the genomic bin coordinates, their read count and state classification.
#' }
#' \item{segments}{
#' A \link{GRanges-class} object containing regions and their state classification.
#' }
#' \item{weights}{Weight for each component.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Estimated parameters of the emission distributions.}
#' \item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#'
#' @seealso findCNVs
NULL

#' Bivariate Hidden Markov Model
#'
#' The \code{aneuBiHMM} object is output of the function \code{\link{findCNVs.strandseq}} and is basically a list with various entries. The class() attribute of this list was set to "aneuBiHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
#'
#' @name aneuBiHMM
#' @return
#' \item{ID}{An identifier that is used in various \pkg{\link{AneuFinder}} functions.}
#' \item{bins}{
#' A \link{GRanges-class} object containing the genomic bin coordinates, their read count and state classification.
#' }
#' \item{segments}{
#' A \link{GRanges-class} object containing regions and their state classification.
#' }
#' \item{weights}{Weight for each component.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Estimated parameters of the emission distributions.}
#' \item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#'
#' @seealso findCNVs.strandseq
NULL


