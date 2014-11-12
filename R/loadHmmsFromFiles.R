loadHmmsFromFiles <- function(hmm.list) {

	## Intercept user input
	if (check.univariate.modellist(hmm.list)!=0) {
		cat("loading univariate HMMs from files ...")
		mlist <- NULL
		for (modelfile in hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		hmm.list <- mlist
		remove(mlist)
		if (check.univariate.modellist(hmm.list)!=0) stop("argument 'hmm.list' expects a list of univariate hmms or a list of files that contain univariate hmms")
		cat(" done\n")
	}
	
	return(hmm.list)

}


