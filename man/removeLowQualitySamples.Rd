\name{removeLowQualitySamples}
\alias{removeLowQualitySamples}

\title{Quality control of SnpSetIllumina objects}
\description{
	Remove samples from a SnpSetIllumina object that show a low quality 
}
\usage{
removeLowQualitySamples(object, min.intensity = 1500, min.gt = 100, subsample = "OPA")
}
\arguments{
  \item{object}{\code{\link{SnpSetIllumina-class}} object}
  \item{min.intensity}{\code{numeric}. Samples that show a median intensity below this 
	     value in either Red or Green channel are removed}
  \item{min.gt}{\code{numeric}. Samples that have less than this amount of valid 
    genotypes are removed}
  \item{subsample}{factor or column name in \code{reporterInfo} slot of \code{object}}
}
\value{
This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}
\examples{
  data(chr17.260)
  chr17.260<-removeLowQualitySamples(chr17.260,min.gt=10)
}
\keyword{manip}
