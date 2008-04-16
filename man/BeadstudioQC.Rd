\name{BeadstudioQC}
\alias{BeadstudioQC}
\alias{pdfBeadstudioQC}
\title{Quality control of Beadstudio report files}
\description{
  When data has been imported using a Beadstudio samplesheet and reportfile,
  these functions can be used to generate quality measures
}
\usage{
BeadstudioQC(object, QClist = list(), arrayType = "Sentrix96")
pdfBeadstudioQC(QClist, basename = "beadstudio", by = 10)
}

\arguments{
  \item{object}{\code{SnpSetIllumina} object.}
  \item{QClist}{list, result of previous call to \code{BeadstudioQC}}
  \item{arrayType}{character, type of array}
  \item{basename}{character, prefix for PDF files. This name will be added before the Barcode of the chip}
  \item{by}{integer, number of samples in barplot, see \code{\link{reportSamplePanelQC}}}
}
\value{
  The \code{BeadstudioQC} function generates a list of \code{\link{QCIllumina}} objects
  The \code{pdfBeadstudioQC} function generates a pdf-file for each \code{QCIllumina} object in the list
}
\author{J. Oosting}
\seealso{\code{\link{pdfQC}},\code{\link{calculateQCarray}}}
\keyword{manip}
\keyword{hplot}
