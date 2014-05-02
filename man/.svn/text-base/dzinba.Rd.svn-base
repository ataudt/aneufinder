\name{dzinba}
\alias{dzinba}
\alias{pzinba}
\alias{rzinba}
\title{
Zero inflated negative binomial distribution
}
\description{
  Density, CDF and random number generation for the zero inflated
  negative binomial distribution.
}
\usage{
dzinba(x, size, mu, beta, log = FALSE)
pzinba(x, size, mu, beta, lower.tail = TRUE)
rzinba(n, size, mu, beta)
}
\arguments{
  \item{x}{
count data
}
  \item{size}{
size parameter of the negative binomial distribution
}
  \item{mu}{
mean parameter of the negative binomial distribution
}
  \item{beta}{
zero inflation parameter
}
  \item{log}{
logical indicating whether to return log density values
}
\item{lower.tail}{
logical indicating whether to P(X<=k) or P(X>k)
}
\item{n}{
number of random samples from the distribution
}
}
\details{
  The density of the zero inflated negative binomial distribution is defined as
  P(x) = I(x == 0) * beta + I(x != 0) * (1 - beta) * dnbinom(x, size,
  mu). 
}
\value{
dzinba is the density function, pzinba is the CDF, rzinba returns
\code{n} random samples from the distribution.
}
\references{
}
\author{
Matthias Heinig
}

\seealso{
  \code{\link{fitzinba}}
}
\examples{

x = rzinba(1000, mu=10, size=15, beta=0.1)
hist(x, breaks=50)
lines(0:max(x), dzinba(0:max(x), mu=10, size=15, beta=0.1) * length(x))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
