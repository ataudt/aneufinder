\name{univariate.from.binned.data}
\alias{univariate.from.binned.data}
\title{
Perform the univariate HMM with binned data
}
\description{
Use previously binned data to estimate a Hidden Markov Model (HMM). Output is a model which contains the posterior probabilities for each bin.
}
\usage{
univariate.from.binned.data(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=2, output.if.not.converged=FALSE)
}
\arguments{
\item{binned.data}{
A data frame with 4 columns: 'chromosome', 'start', 'end', 'reads'. The columns 'chromosome', 'start' and 'end' contain the coordinates for consecutive bins on the genome. Column 'reads' contains the number of aligned reads in that bin.

\tabular{rrrr}{
chromosome \tab start \tab end \tab reads\cr
chr1 \tab 0 \tab 199 \tab 0\cr
chr1 \tab 200 \tab 399 \tab 10\cr
chr1 \tab 400 \tab 599 \tab 9\cr
... \tab ... \tab ... \tab ...\cr
}
}
\item{eps}{
Convergence threshold for the Baum-Welch algorithm.
}
\item{max.time}{
The maximum running time for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. Negative numbers are interpreted as no limit.
}
\item{max.it}{
The maximum number of iterations for the Baum-Welch algorithm. Negative numbers are interpreted as no limit.
}
\item{num.trials}{
The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output.
}
\item{eps.try}{
If option num.trials is set to greater than 1, 'eps.try' is used for the trial runs. If unset, 'eps' is used.
}
\item{num.threads}{
Number of threads to use. The default of 2 threads yields the best performance, given there are at least two cores on your machine.
}
\item{output.if.not.converged}{
If set to 'FALSE', no output will be given in the case that the Baum-Welch did not converge to the desired 'eps'. If set to 'TRUE', a model will be returned even if the Baum-Welch did not converge to the desired 'eps'.
}
}

\details{
The function expects a data.frame with 4 columns (see arguments). Output of the function is a list object which contains the posteriors of the Baum-Welch likelihood maximization and the parameters for the Negative Binomial distributions. It also contains miscellaneous other attributes. All attributes can be accessed with the '$' or '[[ ]]' operator. To see which attributes are available, type "names(output)" (without quotation marks).
}
\value{
The output is a list object with various entries:
\item{coordinates}{A data.frame with genomic coordinates. These are the same as in the input.}
\item{reads}{Number of aligned reads in the bins. Same as in input.}
\item{posteriors}{Posterior probabilities of the Baum-Welch estimation.}
\item{loglik}{Logarithm of the likelihood estimator (goodness of the fit).}
\item{iteration}{Number of iterations performed by the Baum-Welch.}
\item{time.in.sec}{Time in seconds spent in the Baum-Welch.}
\item{delta.loglik}{Change of the log(likelihood) in the last step. This should be smaller than the 'epsilon' in case of convergence.}
\item{epsilon}{Convergence threshold used for the Baum-Welch.}
\item{proba}{Estimated posteriors of the first bin.}
\item{A}{Estimated transition matrix.}
\item{distributions}{Estimated distribution parameters.}
\item{softweights}{Weights of the states. 'Soft' means that only posteriors where used to calculate the weights, as opposed to 'hard' assignments, where each bin is assigned to a state by using a threshold.}
\item{proba.initial}{Initial probabilities of the first bin.}
\item{A.initial}{Initial transition matrix.}
\item{distributions.initial}{Initial distribution parameters.}
\item{num.threads}{Number of threads used for the computations.}
}
\keyword{Hidden Markov Model}
\keyword{Baum-Welch}
\keyword{univariate}
\references{}
\author{Aaron Taudt and Maria Colome Tatche}
\note{}
\seealso{}
\examples{
# Load the example dataset with genome coordinates and number of aligned reads for each bin
data(binned.data.example.1)
# Do the univariate HMM
model = univariate.from.binned.data(binned.data.example.1, eps=1, max.time=60)
# Look at the output
names(model)
# Select the posteriors
head(model$posteriors)
}
