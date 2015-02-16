\name{Zero-inflated NegBinomial}
\alias{Zero-inflated NegBinomial}
\alias{dzinbinom}
\alias{pzinbinom}
\alias{qzinbinom}
\alias{rzinbinom}
\title{The Zero-inflated Negative Binomial Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the zero-inflated negative binomial distribution with parameters
  \code{w}, \code{size} and \code{prob}.
}
\usage{
dzinbinom(x, w, size, prob, mu)
pzinbinom(q, w, size, prob, mu, lower.tail=TRUE)
qzinbinom(p, w, size, prob, mu, lower.tail=TRUE)
rzinbinom(n, w, size, prob, mu)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.  If \code{length(n) > 1}, the length
    is taken to be the number required.}
	\item{w}{weight of the zero-inflation. \code{0 <= w <= 1}.}
  \item{size}{target for number of successful trials, or dispersion
    parameter (the shape parameter of the gamma mixing distribution).
    Must be strictly positive, need not be integer.}
  \item{prob}{probability of success in each trial. \code{0 < prob <= 1}.}
  \item{mu}{alternative parametrization via mean: see \sQuote{Details}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-inflated negative binomial distribution with \code{size} \eqn{= n} and
  \code{prob} \eqn{= p} has density
  \deqn{
    p(x) = w + (1-w) \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x}{
    \Gamma(x+n)/(\Gamma(n) x!) p^n (1-p)^x}
  for \eqn{x = 0}, \eqn{n > 0} and \eqn{0 < p \le 1}.

  \deqn{
    p(x) = (1-w) \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x}{
    \Gamma(x+n)/(\Gamma(n) x!) p^n (1-p)^x}
  for \eqn{x = 1, 2, \ldots}, \eqn{n > 0} and \eqn{0 < p \le 1}.

}
\value{
  \code{dnbinom} gives the density,
  \code{pnbinom} gives the distribution function,
  \code{qnbinom} gives the quantile function, and
  \code{rnbinom} generates random deviates.

  Invalid \code{w}, \code{size} or \code{prob} will result in return value
  \code{NaN}, with a warning.
}
\source{
	All the functions use \code{\link{dnbinom}} or \code{\link{rnbinom}}.
}
\seealso{
  \link{Distributions} for standard distributions, including
  \code{\link{dbinom}} for the binomial, \code{\link{dnbinom}} for the negative binomial, \code{\link{dpois}} for the
  Poisson and \code{\link{dgeom}} for the geometric distribution, which
  is a special case of the negative binomial.
}
\examples{
## Plot a zero-inflated negative binomial distribution
x <- 0:20
y <- dzinbinom(x, w=0.1, size=5, mu=12)
plot(x, y)
}
\keyword{distribution}
