\name{reportGenotypeSegmentation}
\alias{reportGenotypeSegmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot genomic view}
\description{
}
\usage{
reportGenotypeSegmentation(object, plotRaw = TRUE, subsample = NULL, panels = 0, minProbes = 10, maxY = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{class \code{SnpSetIllumina} after segmentation}
  \item{plotRaw}{logical}
  \item{subsample}{factor}
  \item{panels}{number of panels on a page}
  \item{minProbes}{minimum number of probes for a chromosome within a panel}
  \item{maxY}{maximum value on vertical scale within panels}
  \item{\dots}{arguments are forwrded to \code{plot}}
}
\details{
}
\value{
  this function is used for its side effects
}
\author{Jan Oosting}
\keyword{hplot}
