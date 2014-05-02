\name{callRegions}
\alias{callRegions}
\title{
call regions from a posterior matrix
}
\description{
Call regions using a cutoff on the posterior probabilities.
}
\usage{
callRegions(avg.posterior, cutoff, posterior.col = "avg.posterior", expr.col = "avg.expression")
}

\arguments{
  \item{avg.posterior}{
matrix or data.frame with state posterior probabilities
}
  \item{cutoff}{
cutoff on the posterior probability
}
  \item{posterior.col}{
column index or name of the posterior probability
}
  \item{expr.col}{
column index or name of the gene expression column
}
}

\details{
  The function expects a data.frame or matrix containing the position
  of bins in the genome, the results of the HMM analysis and if
  available gene expression data:
  chrom (chromosome),
  start (start coordinate of the bin),
  end (end coordinate of the bin),
  posterior (posterior probability a HMM state)
  expr (expression values assigned to each bin)

  The first 3 columns have to be named exactly as given above, the
  colnames of the posterior prob and expression data are taken from the
  arguments to the function.

  The regions will be called by applying the specified cutoff and
  joining all consecutive bins. Posterior probabilities as well as
  expression values will be avaraged over the called region.

}
\value{
A GRanges object with metadata columns  
\item{avg.posterior}{Average state posterior in the region}
\item{avg.expr}{Average gene expression in the region}

}
\references{
}
\author{
Matthias Heinig
}
\note{
}

\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (avg.posterior, cutoff, posterior.col = "avg.posterior", 
    expr.col = "avg.expression") 
{
    require(GenomicRanges)
    called.regions = avg.posterior[which(avg.posterior[, posterior.col] > 
        cutoff), ]
    called.granges = GRanges(seqnames = called.regions$chrom, 
        ranges = IRanges(start = called.regions$start, end = called.regions$end))
    elementMetadata(called.granges) = DataFrame(called.regions[, 
        c(posterior.col, expr.col)])
    reduced.granges = reduce(called.granges)
    called2reduced = matchMatrix(findOverlaps(called.granges, 
        reduced.granges))
    reduced.posterior = tapply(elementMetadata(called.granges)[, 
        posterior.col], called2reduced[, 2], mean, na.rm = T)
    reduced.expression = tapply(elementMetadata(called.granges)[, 
        expr.col], called2reduced[, 2], mean, na.rm = T)
    elementMetadata(reduced.granges) = DataFrame(avg.posterior = reduced.posterior, 
        avg.expression = reduced.expression)
    return(reduced.granges)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
