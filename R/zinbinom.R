

#' The Zero-inflated Negative Binomial Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the zero-inflated negative binomial distribution with parameters
#' \code{w}, \code{size} and \code{prob}.
#' 
#'  The zero-inflated negative binomial distribution with \code{size} \eqn{= n} and
#'  \code{prob} \eqn{= p} has density
#'  \deqn{
#'    p(x) = w + (1-w) \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x}{
#'    w + (1-w) * \Gamma(x+n)/(\Gamma(n) x!) p^n (1-p)^x}
#'  for \eqn{x = 0}, \eqn{n > 0}, \eqn{0 < p \le 1} and \eqn{0 \le w \le 1}.
#'
#'  \deqn{
#'    p(x) = (1-w) \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x}{
#'    (1-w) * \Gamma(x+n)/(\Gamma(n) x!) p^n (1-p)^x}
#'  for \eqn{x = 1, 2, \ldots}, \eqn{n > 0}, \eqn{0 < p \le 1} and \eqn{0 \le w \le 1}.
#'
#' @name zinbinom
#' @return dzinbinom gives the density, pzinbinom gives the distribution function, qzinbinom gives the quantile function, and rzinbinom generates random deviates.
#' @author Matthias Heinig, Aaron Taudt
#' @seealso   \link{Distributions} for standard distributions, including
#'  \code{\link{dbinom}} for the binomial, \code{\link{dnbinom}} for the negative binomial, \code{\link{dpois}} for the
#'  Poisson and \code{\link{dgeom}} for the geometric distribution, which
#'  is a special case of the negative binomial.
NULL

#' @describeIn zinbinom gives the density
#' @param x Vector of (non-negative integer) quantiles.
#' @param w Weight of the zero-inflation. \code{0 <= w <= 1}.
#' @param size Target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob Probability of success in each trial. \code{0 < prob <= 1}.
#' @param mu Alternative parametrization via mean: see \sQuote{Details}.
#' @importFrom stats dnbinom
dzinbinom = function(x, w, size, prob, mu) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(x)))
	}
	regular.density = stats::dnbinom(x, size=size, prob=prob, mu=mu)
	zi = rep(0, length(x))
	zi[x==0] = w
	density = zi + (1-w)*regular.density
		
	return(density)
}

#' @describeIn zinbinom gives the cumulative distribution function
#' @inheritParams dzinbinom
#' @param q Vector of quantiles.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @importFrom stats stepfun
pzinbinom = function(q, w, size, prob, mu, lower.tail=TRUE) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(q)))
	}
	# Use stepfun to make also non-integer values possible
	x = 0:max(q)
	cdf = stats::stepfun(x, c(0, cumsum(dzinbinom(x, w, size, prob, mu))))
	p = cdf(q)
	p[p>1] = 1

	if (lower.tail==TRUE) {
		return(p)
	} else {
		return(1-p)
	}
}
	
#' @describeIn zinbinom gives the quantile function
#' @inheritParams dzinbinom
#' @inheritParams pzinbinom
#' @param p Vector of probabilities.
#' @importFrom stats stepfun
qzinbinom = function(p, w, size, prob, mu, lower.tail=TRUE) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(p)))
	}
	x = 0:5000
	cdf = pzinbinom(x, w, size, prob, mu)
	inv.cdf = stats::stepfun(cdf[-length(cdf)], x)
	q = inv.cdf(p)

	if (lower.tail==TRUE) {
		return(q)
	} else {
		return(rev(q))
	}
}

#' @describeIn zinbinom random number generation
#' @inheritParams dzinbinom
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @importFrom stats rnbinom
rzinbinom = function(n, w, size, prob, mu) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(n)))
	}
	r = stats::rnbinom(n, size=size, prob=prob, mu=mu)
	if (length(n)==1) {
		i = sample(n, round(w * n))
		if (n < 10) warning("n < 10: Random numbers may not be reliable")
	} else {
		i = sample(n, round(w * length(n)))
		if (length(n) < 10) warning("length(n) < 10: Random numbers may not be reliable")
	}
	r[i] = 0
	
	return(r)
}

