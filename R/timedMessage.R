

startTimedMessage <- function(...) {

	x <- paste0(..., collapse='')
	message(x, appendLF=FALSE)
	ptm <- proc.time()
	return(ptm)

}


stopTimedMessage <- function(ptm) {

	time <- proc.time() - ptm
	message(" ", round(time[3],2), "s")

}
