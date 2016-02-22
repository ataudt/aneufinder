

#' Initialize state factor levels and distributions
#'
#' Initialize the state factor levels and distributions for the specified states.
#'
#' @param states A subset of \code{c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy","+10-somy")}.
#' @return A \code{list} with $labels, $distributions and $multiplicity values for the given states.
initializeStates <- function(states) {

	possible.states <- c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy","5-somy","6-somy","7-somy","8-somy","9-somy","+10-somy")
	possible.distributions <- factor(c("zero-inflation"='delta',
																			"0-somy"='dgeom',
																			"1-somy"='dnbinom',
																			"2-somy"='dnbinom',
																			"3-somy"='dnbinom',
																			"4-somy"='dnbinom',
																			"5-somy"='dnbinom',
																			"6-somy"='dnbinom',
																			"7-somy"='dnbinom',
																			"8-somy"='dnbinom',
																			"9-somy"='dnbinom',
																			"+10-somy"='dnbinom'), levels=c('delta','dgeom','dnbinom','dbinom'))
	possible.multiplicity <- c("zero-inflation"=0,
										"0-somy"=0,
										"1-somy"=1,
										"2-somy"=2,
										"3-somy"=3,
										"4-somy"=4,
										"5-somy"=5,
										"6-somy"=6,
										"7-somy"=7,
										"8-somy"=8,
										"9-somy"=9,
										"+10-somy"=10)
	if (any(!(states %in% possible.states))) {
		stop('argument \'states\' accepts only entries from c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy","5-somy","6-somy","7-somy","8-somy","9-somy","+10-somy")')
	}
	labels <- factor(states, levels=possible.states[possible.states %in% states])
	distributions <- possible.distributions[states]
	multiplicity <- possible.multiplicity[states]
	
	# Return list
	l <- list(labels=labels, distributions=distributions, multiplicity=multiplicity)
	return(l)
}


