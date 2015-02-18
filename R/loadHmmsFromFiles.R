#' Load HMMs from files
#'
#' Load \code{\link{aneuHMM}} objects from file into a list.
#'
#' @param hmm.list A list of files that contain \code{\link{aneuHMM}} objects.
#' @param strict If any of the loaded objects is not a \code{\link{aneuHMM}} object, an error (\code{strict=TRUE}) or a warning (\code{strict=FALSE}) will be generated.
#' @return A list() containing all loaded \code{\link{aneuHMM}} objects.
#' @author Aaron Taudt
#' @export
loadHmmsFromFiles <- function(hmm.list, strict=FALSE) {

	if (is.hmm(hmm.list)) {
		return(hmm.list)
	} else if (is.character(hmm.list)) {
		message("loading univariate HMMs from files ...", appendLF=F); ptm <- proc.time()
		mlist <- list()
		for (modelfile in hmm.list) {
			mlist[[modelfile]] <- get(load(modelfile))
			if (!is.hmm(mlist[[modelfile]])) {
				if (strict) {
					time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
					stop("File ",modelfile," does not contain an ",class.univariate.hmm," object.")
				} else {
					class(modelfile) <- class.univariate.hmm
					warning("File ",modelfile," does not contain an ",class.univariate.hmm," object. Class attribute corrected.")
				}
			}
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		return(mlist)
	} else if (is.list(hmm.list)) {
		index <- which(unlist(lapply(hmm.list, function(hmm) { !is.hmm(hmm) })))
		if (length(index)>0) {
			if (strict) {
				stop("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '))
			} else {
				for (ind in index) {
					class(hmm.list[[ind]]) <- class.univariate.hmm
				}
				warning("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '),". Class attributes corrected.")
			}
		}
		return(hmm.list)
	}
}

is.hmm <- function(hmm) {
	if (class(hmm)==class.univariate.hmm) return(TRUE)
	return(FALSE)
}

