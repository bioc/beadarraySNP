\name{segmentate}
\alias{segmentate}
\title{Segmentation for SnpSetIllumina objects}
\description{
  Use snapCGH package to perform segmentation
}
\usage{
segmentate(object, method = c("DNACopy", "HMM", "BioHMM", "GLAD"), normalizedTo = 2, doLog = TRUE, doMerge = FALSE, subsample = "OPA")
}
\arguments{
  \item{object}{class \code{SnpSetIllumina}}
  \item{method}{char, type of segmentation}
  \item{normalizedTo}{numeric}
  \item{doLog}{logical, perform transformation before segmentation, see \code{convert2seglist}}
  \item{doMerge}{logical, perform merging of close states}
  \item{subsample}{factor}
}
\details{
}
\value{
  SnpSetIllumina object with elements \code{observed, states and predicted} set in the \code{AssayData} slot
}
\author{Jan Oosting}
\keyword{manip}
