\name{Convert to GenomicRanges}
\alias{hmm2GRanges}
\alias{hmmList2GRanges}
\alias{binned2GRanges}

\title{
Convert various \code{\link{aneufinder}} objects into \code{\link{GenomicRanges}} objects.
}
\description{
Convert \code{\link{hmm}}s to \code{\link{GRanges}} objects. Only coordinates and information per position (e.g. reads) are converted, other information specific to the HMM (e.g. loglikelihood, model parameters) are not converted.
}
\usage{
hmm2GRanges(hmm, reduce=TRUE)

hmmList2GRangesList(hmm.list, reduce=TRUE, numCPU=1, consensus=FALSE)

binned2GRanges(binned.data, chrom.length.file=NULL, offset=0)
}
\arguments{
\item{hmm.list}{A list with \code{\link{hmm}}s or a list of files that contain \code{\link{hmm}}s.}
\item{hmm}{An \code{\link{hmm}} object, produced by \code{\link{findCNVs}}.}
\item{reduce}{Either TRUE or FALSE, indicating whether consecutive bins with the same state should be merged into bigger regions.}
\item{numCPU}{Number of CPUs that are used for computations. Requires \option{library(doParallel)}.}
\item{binned.data}{A data.frame with coordinates and binned read counts.}
\item{chrom.length.file}{
A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.
}
\item{offset}{This offset will be added to the coordinates.}
}

\value{
Return value is a \code{\link{GRanges}} object or a \code{\link{GRangesList}} object, respectively.
}
\keyword{GRanges}
\author{Aaron Taudt}
\examples{
}
