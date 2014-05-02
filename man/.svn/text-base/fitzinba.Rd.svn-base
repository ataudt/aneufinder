\name{fitzinba}
\alias{fitzinba}

\title{
Fit parameters of a zero inflated negative binomial distribution
}
\description{
Fit parameters of a zero inflated negative binomial distribution
}
\usage{
fitzinba(x, weight = 1)
}
\arguments{
  \item{x}{vector of count data}
  \item{weight}{weight of each data point (used in the EM algorithm)}
}
\details{
  The density of the zero inflated negative binomial distribution is defined as
  P(x) = I(x == 0) * beta + I(x != 0) * (1 - beta) * dnbinom(x, size,
  mu). The maximum likelihood estimate of the parameters is found using
  numerical optimization (see \code{\link{optim}}).
}
\value{
  Returns a vector with three named values:
  \item{size}{size parameter}
  \item{mu}{mean parameter}
  \item{beta}{zero inflation parameter}
}
\references{
}
\author{
Matthias Heinig
}

\seealso{
\code{\link{dzinba}}, \code{\link{pzinba}}, \code{\link{rzinba}}
}
\examples{
x = rzinba(n=1000, beta=0.1, mu=10, size=15)
fitzinba(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

