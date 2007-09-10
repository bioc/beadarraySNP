\name{plotGoldenGate4OPA}
\alias{plotGoldenGate4OPA}
\title{Plot Golden Gate genomic view}
\description{
  Plots a full genome view based on 4 subsamples of Illumina Golden Gate data
}
\usage{
plotGoldenGate4OPA(object, cn.sum = NULL, sample = 1, plotRaw = FALSE, main = NULL, interact = FALSE, ...)
}
\arguments{
  \item{object}{class \code{SnpSetIllumina}}
  \item{cn.sum}{list containing genomic states, see \code{createCNSummary} }
  \item{sample}{identifier to select the sample within the object}
  \item{plotRaw}{logical, plot raw data points}
  \item{main}{charr, Title of plot}
  \item{interact}{logical, plot should be usable for interactive copy number determination \code{interactiveCNselect}}
  \item{\dots}{extra arguments are formwarded to \code{plot}}
}
\details{
}
\value{
  list, see \code{\link{createCNSummary}}
}
\author{Jan Oosting}
\seealso{\code{\link{segmentate}, \link{alterCN}, \link{interactiveCNselect} \link{createCNSummary}} }
\keyword{hplot}
\keyword{iplot}
