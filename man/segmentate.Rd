\name{segmentate}
\alias{segmentate}
\title{Segmentation for SnpSetIllumina objects}
\description{
  Use snapCGH package to perform segmentation
}
\usage{
segmentate(object, method = c("DNACopy", "HMM", "BioHMM", "GLAD"), normalizedTo = 2, doLog = TRUE,
           doMerge = FALSE, useLair = FALSE, subsample = "OPA", alpha = 0.01)
}
\arguments{
  \item{object}{class \code{SnpSetIllumina}}
  \item{method}{char, type of segmentation}
  \item{normalizedTo}{numeric}
  \item{doLog}{logical, perform transformation before segmentation, see \code{convert2seglist}}
  \item{doMerge}{logical, perform merging of close states}
  \item{useLair}{logical, Also segmentate on lair}
  \item{subsample}{factor}
  \item{alpha}{numeric, probability threshold to distinguish segments}
}
\value{
  SnpSetIllumina object with elements \code{observed, states and predicted} set in the \code{AssayData} slot
}
\author{Jan Oosting}
\keyword{manip}
