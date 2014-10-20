\name{aneufinder HMM object}
\alias{aneufinder.hmm}
\alias{hmm}

\title{HMM object}
\description{
The HMM object is output of the function \code{\link{findCNVs}} and is basically a list with various entries. The class() attribute of this list was set to "aneufinder.hmm". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
}
\value{
\item{reads}{Vector of read counts for each bin.}
\item{num.bins}{Number of bins. Same as \option{length(hmm$reads)}.}
\item{num.states}{Number of states that define this HMM.}
\item{num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \option{eps}.}
\item{time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \option{eps}.}
\item{loglik.delta}{Change in loglikelihood for the last iteration (should be smaller than \option{eps})}
\item{states}{Vector of states for each bin. Type \option{levels(hmm$states)} to get the states that were used.}
\item{A}{Matrix of transition probabilities for each state into each state.}
\item{proba}{Probabilities for the first bin.}
\item{loglik}{Final loglikelihood after the last iteration.}
\item{weights}{Weight for each component.}
\item{ini.proc}{Initialization procedure used.}
\item{A.initial}{Initial \option{A} at the beginning of the Baum-Welch.}
\item{proba.initial}{Initial \option{proba} at the beginning of the Baum-Welch.}
\item{num.threads}{Number of threads that were used in the Baum-Welch.}
\item{error}{This is 0 unless an error occurred during the Baum-Welch.}
\item{read.cutoff}{The read cutoff as specified by the option \option{read.cutoff.quantile}.}
\item{eps}{Convergence threshold for the Baum-Welch.}
\item{distributions}{Estimated parameters of the employed distributions. See \code{\link{findCNVs}} for a description of which distribution is used for which state.}
\item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
\item{ID}{User specified identifier for this object.}
\item{coordinates}{Genomic coordinates of the bins.}
\item{seqlengths}{The lengths of the chromosomes.}
\item{filter.reads}{Indicates if option \option{filter.reads} was used.}
}

\keyword{hmm}
\seealso{\code{\link{hmm2GRanges}}, \code{\link{export.hmms2bed}}}

