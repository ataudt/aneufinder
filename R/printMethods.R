#' Print aneuHMM object
#' 
#' @param x An \code{\link{aneuHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.aneuHMM <- function(x, ...) {
    
    message("$ID")
    print(x$ID)
    message("\n$segments")
    print(x$segments)
    message("\nUse the list operator $ to access all elements of this object.")
  
}


#' Print aneuBiHMM object
#' 
#' @param x An \code{\link{aneuBiHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.aneuBiHMM <- function(x, ...) {
    
    message("$ID")
    print(x$ID)
    message("\n$segments")
    print(x$segments)
    message("\nUse the list operator $ to access all elements of this object.")
  
}

