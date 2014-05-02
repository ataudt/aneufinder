\name{zinba.hmm.posterior}
\alias{zinba.hmm.posterior}
\title{
Compute state posterior probabilities for a HMM with zinba emission
distributions.
}
\description{
Compute state posterior probabilities for a HMM with zinba emission
distributions. First the intial and transition probabilities of the
model are estimated using the Baum-Welch algorithm (considering the
parameters of the emission distributions fixed). Then the state
posteriors are computed using the forward-backward algorithm.
}
\usage{
zinba.hmm.posterior(x, mu, size, beta, max.it = 500, eps = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
count data vector
}
  \item{mu}{
vector of length \code{p} containing the \code{mu} parameters of all
\code{p} states.
}
  \item{size}{
vector of length \code{p} containing the \code{size} parameters of all
\code{p} states.
}
  \item{beta}{
vector of length \code{p} containing the \code{beta} parameters of all
\code{p} states.
}
  \item{max.it}{
maximum number of iterations
}
  \item{eps}{
stop iterations when the gain in loglikelihood is less than \code{eps}
}
}
\details{

}
\value{
  A matrix of posterior probabilities for each data point in the
  sequence (rows) and each state (cols). The matrix also has attributes
  accessible by the \code{attr} function:
  \item{A}{state transition matrix estimated by Baum-Welch}
  \item{proba}{initial state probabilities}
}
\references{
}
\author{
Matthias Heinig
}
\note{
}


\seealso{
  \code{\link{zinbacopula.hmm.posterior}}
  
}
\examples{
# define model parameters for two states
mu = c(5, 10)
size = c(10, 10)
beta = c(0.1, 0.02)

# generate observation sequences
x = c(rzinba(1000, mu=mu[1], size=size[1], beta=beta[1]),
  rzinba(1000, mu=mu[2], size=size[2], beta=beta[2]),
  rzinba(1000, mu=mu[1], size=size[1], beta=beta[1]))

# get the posteriors
posterior = zinba.hmm.posterior(x, mu, size, beta)

# plot the sequence with state assignments
plot(x, col=apply(posterior, 1, which.max))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
