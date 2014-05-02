\name{Zinba-class}
\Rdversion{1.1}
\docType{class}
\alias{Distribution-class}
\alias{LogNormal-class}
\alias{Zinba-class}
\alias{fit,Zinba,numeric,numeric-method}
\alias{logLikelihood,Zinba,numeric-method}

\title{Class \code{Distribution}}
\description{
Classes to represent distributions in mixture models in particular to be
used with the \code{\link{em}} function. Currently Zinba is used but
LogNormal is also available.
}
\section{Objects from the Class Zinba}{
Objects can be created by calls of the form \code{new("Zinba",
  ...)}. For a definition of the density see \code{\link{dzinba}}.
}
\section{Slots}{
  \describe{
    \item{\code{size}:}{size parameter of the negative binomial distribution}
    \item{\code{mu}:}{mean parameter of the negative binomial distribution}
    \item{\code{beta}:}{zero inflation parameter}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Distribution}"}, directly.
}
\section{Methods}{
  \describe{
    \item{fit}{\code{signature(model = "Zinba", x = "numeric",
	logExpectation = "numeric")}: used in the maximization step }
    \item{logLikelihood}{\code{signature(model = "Zinba", x =
	"numeric")}: used in the expectation step }
  }
}

\section{Objects from the Class LogNormal}{
Objects can be created by calls of the form \code{new("LogNormal", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{meanlog}:}{mean parameter (log space)}
    \item{\code{sdlog}:}{standard deviation parameter (log space)}
  }
}


\author{
Matthias Heinig
}

\seealso{
\code{\link{dzinba}}, \code{\link{dlnorm}}
}
\examples{
showClass("Zinba")
}
\keyword{classes}
