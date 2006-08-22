\name{pdfChromosomesSmoothCopyNumber}
\alias{pdfChromosomesSmoothCopyNumber}
\alias{pdfSamplesSmoothCopyNumber}

\title{reportWrappers}
\description{
  Functions that help create pdf reports
}
\usage{
pdfChromosomesSmoothCopyNumber(object, filename, ...)
pdfSamplesSmoothCopyNumber(object, filename, ...)
}

\arguments{
  \item{object}{SnpSetIllumina object}
  \item{filename}{filename of output pdf file}
  \item{\dots}{arguments for report functions}
}
\details{
  These functions set up and perfom reporting to pdf files.
}
\value{
  This function is used for its side effects
}
\author{ Jan Oosting}
\seealso{ \code{\link{reportChromosomesSmoothCopyNumber}}, \code{\link{reportSamplesSmoothCopyNumber}} }
\examples{
  \dontrun{data(chr17.260)}
  \dontrun{data.nrm<-standardNormalization(chr17.260)}
  \dontrun{pdfChromosomesSmoothCopyNumber(data.nrm, "Chr17.pdf", grouping=pData(data.nrm)$Group,smooth.lambda=4)}
}
\keyword{hplot}
