\name{dzinba.copula}
\alias{dzinba.copula}
\alias{pzinba.copula}
\alias{rzinba.copula}
\title{
Density, CDF and random number generation from a bivariate zinba copula
}
\description{
Density, CDF and random number generation from a bivariate zinba copula
distribution
}
\usage{
dzinba.copula(x, y, fit, log = FALSE)
pzinba.copula(x, y, fit)
rzinba.copula(n, fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
count data vector
}
  \item{y}{
count data vector
}
  \item{fit}{
parameters of the zinba copula distribution as returned by fit.zinba.copula:
  marginal.x (named parameter vector as returned by \code{\link{fitzinba}}),
  marginal.y (named parameter vector as returned by \code{\link{fitzinba}}), 
  sigma (covariance matrix of the transformed variables) and 
  mu (mean vector of the transformed variables)
}
  \item{log}{
logical indicating whether log density should be returned
}
  \item{n}{
number of random samples to be generated
}
}
\details{
dzinba.copula is the density function, zinba.copula is the CDF and
rzinba.copula returns n random samples from the distribution.
}
\value{
}
\references{
}
\author{
Matthias Heinig
}
\note{
}

\seealso{
\code{\link{fit.zinba.copula}}, \code{\link{fitzinba}}
}
\examples{

# define a zinba copula distribution
fit =  list(mu=c(0, 0),
  marginal.x=c(size=10, mu=10, beta=0.05),
  marginal.y=c(size=10, mu=20, beta=0.02),
  sigma=matrix(c(1, 0.5, 0.5, 1), ncol=2))

# enumerate a 100x100 square
x = 0:60
y = 0:60

# compute the density
d = t(sapply(x, function(xx) sapply(y, function(yy) dzinba.copula(xx, yy, fit))))

# plot the density
image(x=x, y=y, z=log(d))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
