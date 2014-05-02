\name{histoneHMM-package}
\alias{histoneHMM-package}
\alias{histoneHMM}
\docType{package}
\title{Analysis of histone modification data using a HMM}
\description{
  histoneHMM is a software to analyse ChIP-seq data of histone modifications with broad genomic footprints like H3K27me3. It allows for calling modified regions in single samples as well as for calling differentially modified regions in a comparison of two samples.
}
\details{
\tabular{ll}{
Package: \tab histoneHMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-18\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

The most important functions of the package are
\code{\link{zinba.hmm.posterior}},
\code{\link{zinbacopula.hmm.posterior}},
\code{\link{em}},
\code{\link{fitzinba}} and
\code{\link{fit.zinba.copula}}.

}
\author{
Author: Matthias Heinig and Maria Colome Tatche
Maintainer: Matthias Heinig <heinig@molgen.mpg.de>
}
\references{http://histonehmm.molgen.mpg.de}
\keyword{ package }
\seealso{
}
\examples{

# define model parameters for two states
mu = c(5, 10)
size = c(10, 10)
beta = c(0.1, 0.02)

# generate observation sequences
x = c(rzinba(1000, mu=mu[1], size=size[1], beta=beta[1]),
  rzinba(1000, mu=mu[2], size=size[2], beta=beta[2]),
  rzinba(1000, mu=mu[1], size=size[1], beta=beta[1]))

# get the posteriors
posterior = zinba.hmm.posterior(x, mu, size, beta)

# plot the sequence with state assignments
plot(x, col=apply(posterior, 1, which.max))

}
