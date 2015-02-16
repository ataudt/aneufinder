dzinbinom = function(x, w, size, prob, mu) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(x)))
	}
	regular.density = dnbinom(x, size=size, prob=prob, mu=mu)
	zi = rep(0, length(x))
	zi[x==0] = w
	density = zi + (1-w)*regular.density
		
	return(density)
}

pzinbinom = function(q, w, size, prob, mu, lower.tail=TRUE) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(q)))
	}
	# Use stepfun to make also non-integer values possible
	x = 0:max(q)
	cdf = stepfun(x, c(0, cumsum(dzinbinom(x, w, size, prob, mu))))
	p = cdf(q)
	p[p>1] = 1

	if (lower.tail==TRUE) {
		return(p)
	} else {
		return(1-p)
	}
}
	
qzinbinom = function(p, w, size, prob, mu, lower.tail=TRUE) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(p)))
	}
	x = 0:5000
	cdf = pzinbinom(x, w, size, prob, mu)
	inv.cdf = stepfun(cdf[-length(cdf)], x)
	q = inv.cdf(p)

	if (lower.tail==TRUE) {
		return(q)
	} else {
		return(rev(q))
	}
}

rzinbinom = function(n, w, size, prob, mu) {
	if (w < 0 || w > 1) {
		warning("NaNs returned, w needs to be between 0 and 1")
		return(rep(NaN, length(n)))
	}
	r = rnbinom(n, size=size, prob=prob, mu=mu)
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

