# aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
# Copyright (C) 2015  Aaron Taudt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Initialize state factor levels and distributions
#'
#' Initialize the state factor levels and distributions for the specified states.
#'
#' @param states A subset of \code{c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy")}.
initializeStates <- function(states) {

	possible.states <- c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy")
	possible.distributions <- factor(c("zero-inflation"='delta',
																			"nullsomy"='dgeom',
																			"monosomy"='dnbinom',
																			"disomy"='dnbinom',
																			"trisomy"='dnbinom',
																			"tetrasomy"='dnbinom',
																			"multisomy"='dnbinom'), levels=c('delta','dgeom','dnbinom','dbinom'))
	possible.multiplicity <- c("zero-inflation"=0,
										"nullsomy"=0,
										"monosomy"=1,
										"disomy"=2,
										"trisomy"=3,
										"tetrasomy"=4,
										"multisomy"=5)
	if (any(!(states %in% possible.states))) {
		stop('argument \'states\' accepts only entries from c("zero-inflation","nullsomy","monosomy","disomy","trisomy","tetrasomy","multisomy")')
	}
	labels <- factor(states, levels=possible.states[possible.states %in% states])
	distributions <- possible.distributions[states]
	multiplicity <- possible.multiplicity[states]
	
	# Return list
	l <- list(labels=labels, distributions=distributions, multiplicity=multiplicity)
	return(l)
}


