\name{plotGoldenGate4OPA}
\alias{plotGoldenGate4OPA}
\alias{plotGenomePanels}

\title{Plot Golden Gate genomic view}
\description{
  Plots a full genome view based on 4 subsamples of Illumina Golden Gate data
}
\usage{
plotGoldenGate4OPA(object, cn.sum = NULL, sample = 1, plotRaw = FALSE, main = NULL, interact = FALSE, ...)
plotGenomePanels(object, cn.sum = NULL, sample = 1, plotRaw = FALSE, main = NULL, interact = FALSE, panels = list(p1 = c(1, 2, 3, 22), p2 = c(5, 6, 7, 8, 9), p3 = c(10, 11, 12, 13, 14, 15, 21), p4 = c(4, 16, 17, 18, 19, 20, 98, 99)), ...)
}
\arguments{
  \item{object}{class \code{SnpSetIllumina}}
  \item{cn.sum}{list containing genomic states, see \code{createCNSummary} }
  \item{sample}{identifier to select the sample within the object}
  \item{plotRaw}{logical, plot raw data points}
  \item{main}{character, Title of plot}
  \item{interact}{logical, plot should be usable for interactive copy number determination \code{interactiveCNselect}}
  \item{panels}{list, vectors of chromosomes for each panel}
  \item{\dots}{extra arguments are formwarded to \code{plot}}
}
\details{
  prepare interactive selection
}
\value{
  list, see \code{\link{createCNSummary}}
}
\author{Jan Oosting}
\seealso{\code{\link{segmentate}, \link{alterCN}, \link{interactiveCNselect} \link{createCNSummary}} }
\keyword{hplot}
\keyword{iplot}
