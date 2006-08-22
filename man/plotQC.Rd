\name{plotQC}
\alias{plotQC}
\title{ Spatial plots of array QC information }
\description{
  Plots array wide summary information using the layout of the physical medium
}
\usage{
plotQC(object, type)
}
\arguments{
  \item{object}{ object that contains QC information. e.g. \code{QCIllumina-class} }
  \item{type}{\code{character}, the type of information to plot, currently the 
    following types are supported: \code{"intensityMed","greenMed","redMed",
    "validn","annotation" and "samples"} }
}
\value{
  The function is used for its side effects
}
\author{Jan Oosting}
\seealso{\code{\link{pdfQC}}, \code{\link{reportSamplePanelQC}}}
\examples{
  data(QC.260)
  plotQC(QC.260,"greenMed")
}
\keyword{hplot}
