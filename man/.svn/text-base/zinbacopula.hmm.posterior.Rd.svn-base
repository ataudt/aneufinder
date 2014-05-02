\name{zinbacopula.hmm.posterior}
\alias{zinbacopula.hmm.posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
State posterior probabilities for a HMM with zinba copula emission probabilities
}
\description{
State posterior probabilities for a HMM with zinba copula emission
probabilities.  First the intial and transition probabilities of the
model are estimated using the Baum-Welch algorithm (considering the
parameters of the emission distributions fixed). Then the state
posteriors are computed using the forward-backward algorithm.
}
\usage{
zinbacopula.hmm.posterior(x, fits, max.it = 500, eps = 1e-04)
}
\arguments{
  \item{x}{
a two row matrix with bivariate observations in each column
}
  \item{fits}{
a list of zinba copula fits as returned by \code{\link{fit.zinba.copula}}
}
  \item{max.it}{
maximum number of iterations
}
  \item{eps}{
stop iterations when the gain in loglikelihood is less than eps
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
\code{\link{zinba.hmm.posterior}}
}
\examples{

# define model parameters for two states
fits = list(state1=list(mu=c(0, 0),
    marginal.x=c(size=10, mu=10, beta=0.05),
    marginal.y=c(size=10, mu=13, beta=0.02),
    sigma=matrix(c(1, 0.5, 0.5, 1), ncol=2)),
  state2=list(mu=c(0, 0),
    marginal.x=c(size=10, mu=13, beta=0.02),
    marginal.y=c(size=10, mu=10, beta=0.05),
    sigma=matrix(c(1, 0.5, 0.5, 1), ncol=2)))


# generate observation sequences
x = rbind(rzinba.copula(1000, fits[[1]]),
  rzinba.copula(1000, fits[[2]]),
  rzinba.copula(1000, fits[[1]]))

true.state = rep(c(1, 2, 1), each=1000)

# get the posteriors
posterior = zinbacopula.hmm.posterior(t(x), fits)

# check the classification
predicted.state = apply(posterior, 1, which.max)
confusion = table(true.state, predicted.state)

# plot the observations and assign colors according to correct classifications
plot(x, col=factor(interaction(true.state, predicted.state)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
