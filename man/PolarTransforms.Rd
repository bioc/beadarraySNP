\name{PolarTransforms}
\alias{RG2polar}
\alias{polar2RG}
\concept{SnpSetIllumina}

\title{Polar transformations}

\description{
   Perform polar transforms on Illumina Golden Gate bead arrays
}

\usage{
   RG2polar(object,trig=FALSE)
   polar2RG(object,trig=FALSE)
}
\arguments{
   \item{object}{SnpSetIllumina object}
   \item{trig}{Logical, use geometric distance intensity. Otherwise use addition of intensities}
}

\details{
  \code{RG2polar} transforms the \code{R} and \code{G} matrices to \code{theta} and
  \code{intensity} matrices. Note that the intensity value is the sum of \code{R}
  and \code{G} and not the geometric distance to the origin.
  
  \code{polar2RG} performs the reverse transformation
}
\value{
  This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}
}

\examples{
  data(chr17.260)
  data.polar<-RG2polar(chr17.260)
  plot(assayData(data.polar)$theta,assayData(data.polar)$intensity)
}
\keyword{manip}
