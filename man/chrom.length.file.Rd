\name{File with chromosome lengths}
\alias{chrom.length.file}

\title{Format of the chrom.length.file}
\description{
The chrom.length.file can or must be used in the following \code{\link{aneufinder}} functions: \code{\link{bed2binned}}, \code{\link{bedGraph2binned}}
}
\section{Format}{
In the first column, it contains the chromosome names and in the second column, it contains the chromosome lengths in basepairs. The file has no header. Example:
\tabular{ll}{
chr1 \tab 249250621 \cr
chr2 \tab 243199373 \cr
chr3 \tab 198022430 \cr
chr4 \tab 191154276 \cr
... \tab ... \cr
}
}
\keyword{chromosome}
\keyword{length}
