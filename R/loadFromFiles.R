#' Load \pkg{AneuFinder} objects from file
#'
#' Wrapper to load \pkg{\link{AneuFinder}} objects from file and check the class of the loaded objects.
#'
#' @param files A list of \code{\link{GRanges}}, \code{\link{GRangesList}}, \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects or a character vector with files that contain such objects.
#' @param check.class Any combination of \code{c('GRanges', 'GRangesList', 'aneuHMM', 'aneuBiHMM')}. If any of the loaded objects does not belong to the specified class, an error is thrown.
#' @return A list of \code{\link{GRanges}}, \code{\link{GRangesList}}, \code{\link{aneuHMM}} or \code{\link{aneuBiHMM}} objects.
#' @export
#' @examples
#'## Get some files that you want to load
#'folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'files <- list.files(folder, full.names=TRUE)
#'## Load and plot the first ten
#'hmms <- loadFromFiles(files[1:10])
#'lapply(hmms, plot, type='profile')
#'
loadFromFiles <- function(files, check.class=c('GRanges', 'GRangesList', 'aneuHMM', 'aneuBiHMM')) {

    # ptm <- startTimedMessage("Loading data from files ...")
    if (is.null(files)) {
        # stopTimedMessage(ptm)
        return(files)
    }
    if (any(! check.class %in% c('GRanges', 'GRangesList', class.univariate.hmm, class.bivariate.hmm))) {
        stop("Argument 'check.class' must contain any combination of c('", paste0(c('GRanges', 'GRangesList', class.univariate.hmm, class.bivariate.hmm), collapse="', '"), "').")
    }
    modellist <- list()
    if (is.character(files)) {
        for (file in files) {
            temp.env <- new.env()
            model <- get(load(file, envir=temp.env), envir=temp.env)
            if (! class(model) %in% check.class) {
                stop("File '", file, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[file]] <- model
        }
    } else if (class(files) %in% check.class) {
        modellist[[1]] <- files
    } else if (is.list(files)) {
        for (file in files) {
            model <- file
            if (! class(model) %in% check.class) {
                stop("List entry '", length(modellist)+1, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[length(modellist)+1]] <- model
        }
        names(modellist) <- names(files)
    } else if (! class(files) %in% check.class) {
        stop("Input does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
    }
    # stopTimedMessage(ptm)
    return(modellist)

}
