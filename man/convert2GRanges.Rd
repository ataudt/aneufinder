\name{Convert to GenomicRanges}
\alias{hmm2GRanges}
\alias{binned2GRanges}

\title{
Convert various \code{\link{aneufinder}} objects into \code{\link{GenomicRanges}} objects.
}
\description{
Convert \code{\link{binned.data}} and \code{\link{uni.hmm}}s to GRanges objects. Only coordinates and information per position (e.g. reads, posteriors) are converted, other information specific to the HMM (e.g. loglikelihood, model parameters) are not converted.
}
\usage{
binned2GRanges(binned.data)

hmm2GRanges(hmm, reduce=TRUE)
}
\arguments{
\item{binned.data}{A data.frame with coordinates and binned read counts (see \code{\link{binned.data}}).}
\item{hmm}{A \code{\link{uni.hmm}} object, produced by \code{\link{find.aneuploidies}}.}
\item{reduce}{Either TRUE or FALSE, indicating whether consecutive bins with the same state should be merged into bigger regions.}
}

\value{
A GRanges object is returned.
}
\keyword{GRanges}
\author{Aaron Taudt}
\examples{
}
