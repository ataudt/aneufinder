

#' Initialize state factor levels and distributions
#'
#' Initialize the state factor levels and distributions for the specified states.
#'
#' @param states A subset of \code{c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy",...)}.
#' @return A \code{list} with $labels, $distributions and $multiplicity values for the given states.
initializeStates <- function(states) {

	somy.states <- grep('somy', states, value=TRUE)
	somy.numbers <- as.integer(sapply(strsplit(somy.states, '-somy'), '[[', 1))
	names(somy.numbers) <- somy.states

	if ("zero-inflation" %in% states) {
    	multiplicity <- c("zero-inflation"=0, somy.numbers)
	} else {
    	multiplicity <- somy.numbers
	}

	levels.distributions <- c('delta','dgeom','dnbinom','dbinom')
	distributions <- rep(NA, length(states))
	names(distributions) <- states
	distributions[states=='zero-inflation'] <- 'delta'
	distributions[states=='0-somy'] <- 'dgeom'
	distributions[(states != 'zero-inflation') & (states != '0-somy')] <- 'dnbinom'

	# if (any(diff(somy.numbers) > 1)) {
	# 	warning("Copy numbers are not consecutive: ", paste0(somy.states, collapse=', '))
	# }
	# Return list
	states <- factor(states, levels=states)
	distributions <- factor(distributions, levels=levels.distributions)
	l <- list(states=states, distributions=distributions, multiplicity=multiplicity)
	return(l)
}


