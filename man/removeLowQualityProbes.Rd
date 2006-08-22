\name{removeLowQualityProbes}
\alias{removeLowQualityProbes}
\title{Quality control of SnpSetIllumina objects}
\description{
	Remove probes form a SnpSetIllumina object that show a low quality throughout 
	the experiment
}
\usage{
removeLowQualityProbes(object, cutoff = 0.25)
}

\arguments{
  \item{object}{SnpSetIllumina object}
  \item{cutoff}{numeric}
}
\details{
  Probes that have a median value below \code{cutoff} * median value for the
  whole experiment are deleted from the object.
}
\value{
  SnpSetIllumina object 
}
\author{Jan Oosting}
\keyword{manip}
