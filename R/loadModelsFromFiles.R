

#' Load HMMs from files
#'
#' Load \code{\link{aneuHMM}} objects from file into a list.
#'
#' @param hmms A list of files that contain \code{\link{aneuHMM}} objects.
#' @param strict If any of the loaded objects is not a \code{\link{aneuHMM}} object, an error (\code{strict=TRUE}) or a warning (\code{strict=FALSE}) will be generated.
#' @return A list() containing all loaded \code{\link{aneuHMM}} objects.
#' @author Aaron Taudt
#' @export
loadHmmsFromFiles <- function(hmms, strict=FALSE) {

	if (is.hmm(hmms) | is.bihmm(hmms)) {
		return(list(hmms))
	} else if (is.character(hmms)) {
		message("Loading univariate HMMs from files ...", appendLF=FALSE); ptm <- proc.time()
		mlist <- list()
		for (modelfile in hmms) {
			tC <- tryCatch({
				mlist[[modelfile]] <- get(load(modelfile))
			}, error = function(err) {
				stop(modelfile,'\n',err)
			})
			if (!is.hmm(mlist[[modelfile]]) & !is.bihmm(mlist[[modelfile]])) {
				if (strict) {
					time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
					stop("File ",modelfile," does not contain an ",class.univariate.hmm," object.")
				} else {
					class(mlist[[modelfile]]) <- class.univariate.hmm
					warning("File ",modelfile," does not contain an ",class.univariate.hmm," object. Class attribute corrected.")
				}
			}
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		return(mlist)
	} else if (is.list(hmms)) {
		index <- which(unlist(lapply(hmms, function(hmm) { !is.hmm(hmm) & !is.bihmm(hmm) })))
		if (length(index)>0) {
			if (strict) {
				stop("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '))
			} else {
				for (ind in index) {
					class(hmms[[ind]]) <- class.univariate.hmm
				}
				warning("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '),". Class attributes corrected.")
			}
		}
		return(hmms)
	} else if (is.null(hmms)) {
		return(hmms)
	} else {
		warning("Loaded object is not an HMM.")
		return(hmms)
	}
}

is.hmm <- function(hmm) {
	if (class(hmm)==class.univariate.hmm) return(TRUE)
	return(FALSE)
}

is.bihmm <- function(hmm) {
	if (class(hmm)==class.bivariate.hmm) return(TRUE)
	return(FALSE)
}


#' Load GRanges from files
#'
#' Load \code{\link{GRanges}} objects from file into a list.
#'
#' @param files A list of files that contain \code{\link{GRanges}} objects.
#' @return A list() containing all loaded \code{\link{GRanges}} objects.
#' @author Aaron Taudt
#' @export
loadGRangesFromFiles <- function(files) {

	gr.list <- files
	if (class(gr.list)=='GRanges') {
		return(list(gr.list))
	} else if (is.character(gr.list)) {
		message("Loading GRanges from files ...", appendLF=FALSE); ptm <- proc.time()
		mlist <- list()
		for (modelfile in gr.list) {
			tC <- tryCatch({
				mlist[[modelfile]] <- get(load(modelfile))
			}, error = function(err) {
				stop(modelfile,'\n',err)
			})
			if (class(mlist[[modelfile]])!='GRanges') {
				time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
				stop("File ",modelfile," does not contain a GRanges object.")
			}
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		return(mlist)
	} else if (is.list(gr.list)) {
		index <- which(unlist(lapply(gr.list, function(hmm) { class(hmm)!='GRanges' })))
		if (length(index)>0) {
			stop("The following list entries do not contain GRanges objects: ", paste(index, collapse=' '))
		}
		return(gr.list)
	}
}
