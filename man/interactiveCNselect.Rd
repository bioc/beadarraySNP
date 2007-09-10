\name{interactiveCNselect}
\alias{interactiveCNselect}
\title{ Interactive assignment of copynumbers to genomic segments}
\description{
  This function plots the genomic view of a sample, and allows the assignment of a
  discrete copy number to each segment
}
\usage{
interactiveCNselect(object, sample = 1, dnaIndex)
}
\arguments{
  \item{object}{class \code{SnpSetIllumina} after segmentation}
  \item{sample}{Sample identifier within \code{object}}
  \item{dnaIndex}{numeric, measured DNA index of the sample (1=normal)}
}
\details{
}
\value{
  list, see \code{\link{createCNSummary}}
}
\author{Jan Oosting}
\seealso{ \code{\link{segmentate}, \link{alterCN}, \link{plotGoldenGate4OPA} \link{createCNSummary}} }
\keyword{manip}
\keyword{iplot}
