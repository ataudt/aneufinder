karyotypeMeasures <- function(hmms, per.chrom=TRUE) {

	df <- heatmapAneuploidies(hmms, cluster=FALSE, as.data.frame=TRUE)[,-1]
	S <- nrow(df)
	if (per.chrom) {
		D <- sqrt(apply((df - 2)^2, 2, mean))
		V <- apply(df, 2, var)
		tabs <- apply(df, 2, table)
		H <- S*log(S) - S + unlist(lapply(tabs, function(x) { sum(x-x*log(x)) }))
	} else {
		D <- sqrt(mean((df - 2)^2))
		V <- var(df)
		tab <- table(as.vector(df))
		H <- S*log(S) - S + sum(tab-tab*log(tab))
	}
	karyotypeMeasures <- data.frame(divergenceFromDisomic=D, entropicHeterogeneity=H)
	return(karyotypeMeasures)

}


