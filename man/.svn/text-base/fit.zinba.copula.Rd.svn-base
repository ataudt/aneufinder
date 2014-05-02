\name{fit.zinba.copula}
\alias{fit.zinba.copula}
\title{
Maximum likelihood fit of the zinba copula distribution
}
\description{
Maximum likelihood fit of the 2 dimensional zinba copula distribution
}
\usage{
fit.zinba.copula(x, y, weight = 1, marginal.x = NULL, marginal.y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
count data vector
}
  \item{y}{
count data vector
}
  \item{weight}{
weight vector (can be used in the EM algorithm)
}
  \item{marginal.x}{
parameters of the marginal zinba distribution for \code{x}
}
  \item{marginal.y}{
parameters of the marginal zinba distribution for \code{y}
}
}
\details{
if \code{marginal.x} or \code{marginal.y} are NULL they will be
re-estimated, if they are given only the missing parameters will be
re-estimated. The format of the marginal parameters is explained below.
}
\value{
  a list with the parameters of the two marginal distributions:
  \item{marginal.x}{named parameter vector as returned by
    \code{\link{fitzinba}}}
  \item{marginal.y}{named parameter vector as returned by
    \code{\link{fitzinba}}}
  \item{sigma}{covariance matrix of the transformed variables}
  \item{mu}{mean vector of the transformed variables}
}
\references{
}
\author{
Matthias Heinig
}
\note{
}

\seealso{
\code{\link{fitzinba}}, \code{\link{dzinba.copula}}
}
\examples{
x = rzinba.copula(1000, list(mu=c(0, 0), marginal.x=c(size=10, mu=10, beta=0.1), marginal.y=c(size=10, mu=20, beta=0.02), sigma=matrix(c(1, 0.5, 0.5, 1), ncol=2)))

fit = fit.zinba.copula(x[,1], x[,2])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
