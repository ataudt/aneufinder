\name{findCNVs}
\alias{findCNVs}
\title{
Determine Copy Number Variations (CNVs) for single-cell sequencing data
}
\description{
This function takes a GRanges with binned read counts as input and determines the copy number of each bin.
}
\usage{
findCNVs(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE)
}
\arguments{
\item{binned.data}{A GRanges object containing genomic coordinates and read counts (see \code{\link{binned.data}} for a more detailed description).}
\item{ID}{An ID which will be used to identify the returned object in various analysis and plotting functions.}
\item{eps}{Convergence threshold for the Baum-Welch algorithm.}
\item{init}{
Initialization procedure for the distributions. Currently only 'standard' is implemented. The mean and variance of the data will be used to initialize the disomic distribution.
}
\item{max.time}{The maximum running time for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.}
\item{max.iter}{The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.}
\item{num.trials}{The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output. THIS OPTION IS USELESS AT THE MOMENT, starting values are not random.}
\item{eps.try}{If option num.trials is set to greater than 1, \option{eps.try} is used for the trial runs. If unset, \option{eps} is used.}
\item{num.threads}{Number of threads to use. CURRENTLY NOT USED.}
\item{output.if.not.converged}{
If set to FALSE, no output will be given in the case that the Baum-Welch did not converge to the desired \option{eps}. If set to TRUE, a model will be returned even if the Baum-Welch did not converge to the desired \option{eps}.
}
\item{read.cutoff.quantile}{
Read counts which are above this quantile are replaced by the corresponding value. This filtering can increase the speed of the Baum-Welch and is unlikely to affect the results. It may also be helpful to specify a lower cutoff if you are getting the error message 'A nan occurred during the Baum-Welch! ...'.
}
}

\details{
The Hidden Markov Model which is used to classify the bins uses 7 states: state 'nullsomy' with a delta function as emission densitiy (only zero read counts), 'null-mixed' with a geometric distribution (see \code{\link{dgeom}}), 'monosomy','disomy','trisomy','tetrasomy' and 'multisomy' with negative binomials (see \code{\link{dnbinom}}) as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method.
TODO: insert paper
}
\value{
The output is a \code{\link{aneufinder.univariate.hmm}} list object. All attributes can be accessed with the '$' or '[[ ]]' operator. To see which attributes are available, type "names(output)" (without quotation marks).
}
\keyword{Hidden Markov Model}
\keyword{Baum-Welch}
\keyword{univariate}
\references{}
\author{Aaron Taudt and Maria Colome Tatche}
\note{}
\seealso{}
\examples{
}
