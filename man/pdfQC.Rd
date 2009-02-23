\name{pdfQC}
\alias{pdfQC}

\title{ QCreport }
\description{
  Create PDF file with experimental quality control plots
}
\usage{
pdfQC(object, filename = "arrayQC.pdf", by = 10)
}
\arguments{
  \item{object}{QCIllumina object, or list of QCIllumina objects}
  \item{filename}{character, output pdf filename}
  \item{by}{number of samples in barplot, see \code{\link{reportSamplePanelQC}}}
}
\details{
  This function creates a pdf file with QC information. The first page contains 
  8 \code{plotQC} panels showing the spatial distribution of intensities on a SAM plate.
  The following page(s) contain the output of \code{reportSamplePanelQC}
}
\value{
  A PDF file is produced
}
\author{ Jan Oosting }
\seealso{ \code{\link{plotQC}}, \code{\link{reportSamplePanelQC}}, \code{\link{QCIllumina-class}} }
\keyword{hplot}
