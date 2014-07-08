\name{collapse.bins}
\alias{collapse.bins}
\title{
Collapse consecutive bins
}
\description{
The function will collapse consecutive bins which have, for example, the same combinatorial state.
}
\usage{
collapse.bins(data, column2collapseBy=NULL, columns2sumUp=NULL)
}
\arguments{
\item{data}{
A data.frame containing the genomic coordinates in the first three columns.
}
\item{column2collapseBy}{
The number of the column which will be used to collapse all other inputs. If a set of consecutive bins has the same value in this column, they will be aggregated into one bin with adjusted genomic coordinates.
}
\item{columns2sumUp}{
Numbers of columns that will be summed during the aggregation process.
}
}

\details{
The function is used to decrease the file size for export of results into the BED format. It can also be used to cluster annotation data. The following tables illustrate the principle:

Input data:
\tabular{rrrrrr}{
chrom \tab start \tab end \tab column2collapseBy \tab moreColumns \tab columns2sumUp \cr
chr1  \tab     0 \tab 199 \tab                 2 \tab        1 10 \tab           1 3 \cr
chr1  \tab   200 \tab 399 \tab                 2 \tab        2 11 \tab           0 3 \cr
chr1  \tab   400 \tab 599 \tab                 2 \tab        3 12 \tab           1 3 \cr
chr1  \tab   600 \tab 799 \tab                 1 \tab        4 13 \tab           0 3 \cr
chr1  \tab   800 \tab 999 \tab                 1 \tab        5 14 \tab           1 3 \cr
}
Output data:
\tabular{rrrrrr}{
chrom \tab start \tab end \tab column2collapseBy \tab moreColumns \tab columns2sumUp \cr
chr1  \tab     0 \tab 599 \tab                 2 \tab        1 10 \tab           2 9 \cr
chr1  \tab   600 \tab 999 \tab                 1 \tab        4 13 \tab           1 6 \cr
}

}

\value{
An aggregated data.frame with the same format as the input data.frame will be given as output.
}
\keyword{clustering}
\keyword{collapsing}
\keyword{aggregating}
\references{}
\author{Aaron Taudt}
\note{}
\seealso{}
\examples{
}
