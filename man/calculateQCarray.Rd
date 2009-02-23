\name{calculateQCarray}
\alias{calculateQCarray}
\title{ Retrieve QC information from a SnpSetIllumina object}
\description{
  Retrieves QC and identifying information of Illumina Sentrix arrays. 
}
\usage{
calculateQCarray(object, QCobject = NULL, arrayType="Sentrix96")
}

\arguments{
  \item{object}{\code{SnpSetIllumina} object. Should contain information 
	      of a single Sentrix array and a single type of OPA panel}
  \item{QCobject}{\code{QCIllumina-class} object: If set the information in the object 
	      is amended with data from the SnpSetIllumina object}
  \item{arrayType}{character, see \code{\link{arrayType}}}
}
\details{
  Sample summary values are mapped to the physical layout of the 
  Sentrix array using the \code{Row} and \code{Col} columns of the \code{phenoData} 
  slot. These will be available when \code{\link{read.SnpSetIllumina}} is used 
  to create \code{SnpSetIllumina} objects.\cr
  Use successive calls to \code{calculateQCarray} to process Sentrix arrays 
  with multiple probe panels.
  
  If data is read using a samplesheet that defines manifest files it is possible to handle 
  data with multiple manifests and/or multiple Sentrix arrays
  
}
\value{
  A \code{QCIllumina} object, when multiple arrays were combined a list of \code{QCIllumina} objects 
}
\author{Jan Oosting}
\seealso{
  \code{link{QCIllumina-class}}, \code{link{reportSamplePanelQC}}, \code{link{plotQC}}
}
\examples{
  \dontrun{QC<-calculateQCarray(data.raw1)}
  \dontrun{QC<-calculateQCarray(data.raw2,QC)}
}
\keyword{manip}
