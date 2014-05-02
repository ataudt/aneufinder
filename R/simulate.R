simulate.univariate = function(coordinates, transition, emission, initial=1) {

	# Calculate some variables
	numstates = ncol(transition)
	if (numstates!=3) {
		stop("The transition matrix is expected to have 3 columns and 3 rows")
	}
	numbins = nrow(coordinates)

	## Make state vector from transition matrix
	# Generate sample of random numbers
	rsample = runif(numbins,0,1)
	# Integrate transition matrix by row and add -1 in front
	cumtransition = cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
	# Generate the state vector by going through each state
	cat("Generating states from transition matrix...")
	states = matrix(rep(NA,numstates*numbins), ncol=numstates)
	for (irow in 1:numstates) {
		states[,irow] = findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
	}
	statevec = rep(NA,numbins)
	statevec[1] = initial
	for (ibin in 2:numbins) {
		statevec[ibin] = states[ibin,statevec[ibin-1]]
	}
	cat(" done\n")

	## Make reads from state vector and distributions
	# Generate the read counts by drawing from distribution
	cat("Generating reads from emission parameters and states...")
	reads = rep(NA, numbins)
	numbins.in.state = aggregate(rep(1,length(statevec)), list(state=statevec), sum)
	reads[statevec==1] = 0
	if (!is.na(numbins.in.state[2,'x'])) {
		reads[statevec==2] = rnbinom(numbins.in.state[2,'x'], size=emission[2,'r'], prob=emission[2,'p'])
	}
	if (!is.na(numbins.in.state[3,'x'])) {
		reads[statevec==3] = rnbinom(numbins.in.state[3,'x'], size=emission[3,'r'], prob=emission[3,'p'])
	}
	cat(" done\n")

	# Return the output
	out = list(coordinates = coordinates,
				states = statevec,
				reads = reads,
				transition = transition,
				emission = emission
				)
	return(out)

}



simulate.multivariate = function(coordinates, transition, emissions, weights, sigma, use.states, initial=1) {

	lib = require(mvtnorm)
	if (lib == FALSE) {
		install.packages("mvtnorm")
		library(mvtnorm)
	}

	# Calculate some variables
	numstates = ncol(transition)
	numbins = nrow(coordinates)
	nummod = length(emissions)

	## Make state vector from transition matrix
	# Generate sample of random numbers
	rsample = runif(numbins,0,1)
	# Integrate transition matrix by row and add -1 in front
	cumtransition = cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
	# Generate the state vector by going through each state
	cat("Generating states from transition matrix...")
	states = matrix(rep(NA,numstates*numbins), ncol=numstates)
	for (irow in 1:numstates) {
		states[,irow] = findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
	}
	statevec = rep(NA,numbins)
	statevec[1] = initial
	for (ibin in 2:numbins) {
		statevec[ibin] = states[ibin,statevec[ibin-1]]
	}
	# Replace the states by combinatorial states
	statevec = use.states[statevec]
	cat(" done\n")

	## Make reads from state vector and emission distributions
	reads = matrix(rep(NA, nummod*numbins), ncol=nummod)

	rs = unlist(lapply(emissions, "[", 2:3, 'r'))
	ps = unlist(lapply(emissions, "[", 2:3, 'p'))
	ws = unlist(lapply(weights, "[", 1))
	for (istate in use.states) {
		cat("Generating reads for state ",istate,"\n")
		i = which(use.states==istate)
		
		# Convert istate to binary representation
		binary_state = rev(as.integer(intToBits(istate))[1:nummod])
		binary_stateTF = unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])

		# Draw from the multivariate normal
		n = length(which(statevec==istate))
		if (n == 0) next
		cat("drawing from multivariate normal...             \r")
		z = matrix(rmvnorm(n, mean=rep(0,nummod), sigma=sigma[,,i]), ncol=nummod)
		# Transform to uniform space
		cat("transforming to uniform space...                \r")
		u = matrix(apply(z, 2, pnorm), ncol=nummod)
		# Transform to count space using marginals
		cat("transform to count space...                     \r")
		irs = rs[binary_stateTF]
		ips = ps[binary_stateTF]
		for (imod in 1:nummod) {
			mask = statevec==istate
			if (binary_state[imod] == 0) {
				reads[mask,imod] = qzinbinom(u[,imod], w=ws[imod], size=irs[imod], prob=ips[imod])
			} else {
				reads[mask,imod] = qnbinom(u[,imod], size=irs[imod], prob=ips[imod])
			}
		}
		cat("                                                \r")
			
	}

	# Return the output
	out = list(coordinates = coordinates,
				states = statevec,
				reads = reads,
				emissions = emissions,
				transition = transition,
				sigma = sigma,
				use.states = use.states,
				weights = weights
				)
	return(out)

}
