#' Measures for Karyotype Heterogeneity
#' 
#' Computes measures for karyotype heterogeneity. See the Details section for how these measures are defined.
#'
#' We define \eqn{x} as the vector of copy number states for each chromosome and HMM. The number of HMMs is \eqn{S}.
#' \describe{
#' \item{Divergence from disomic:}{\eqn{D = sqrt(mean((x-2)^2))}}
#' \item{Entropic karyotype heterogeneity:}{\eqn{H = S*log(S) - S + sum(table(x)-table(x)*log(table(x)))}}
#' }
#'
#' @param hmms A list with \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @author Aaron Taudt
#' @export
karyotypeMeasures <- function(hmms) {

	karyoMeasures <- list()
	df <- heatmapAneuploidies(hmms, cluster=FALSE, as.data.frame=TRUE)[,-1]
	## Per chromosome
		D <- sqrt(apply((df - 2)^2, 2, mean))
		V <- apply(df, 2, var)
		tabs <- apply(df, 2, table)
		S <- nrow(df)
		H <- S*log(S) - S + unlist(lapply(tabs, function(x) { sum(x-x*log(x)) }))
		karyoMeasures[['per.chromosome']] <- data.frame(divergenceFromDisomic=D, entropicHeterogeneity=H)
# 	## Whole genome
# 		D <- sqrt(mean((df - 2)^2))
# 		V <- var(df)
# 		tab <- table(unlist(df))
# 		S <- prod(dim(df))
# 		H <- S*log(S) - S + sum(tab-tab*log(tab))
# 		karyoMeasures[['whole.genome']] <- data.frame(divergenceFromDisomic=D, entropicHeterogeneity=H)

	return(karyoMeasures)

}


