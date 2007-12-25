\name{QCaccessors}
\alias{arrayType}
\alias{arrayType<-}
\alias{arrayID}
\alias{arrayID<-}

\title{Accessor methods for QC objects}
\description{
	These generic functions set and retrieve properties of quality control 
	objects like \code{\link{QCIllumina-class}}

}
\usage{
arrayType(object)
arrayType(object)<- value
arrayID(object)
arrayID(object)<- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object, possibly derived from class \code{\link{QCIllumina-class}}.}
  \item{value}{}
}
\details{
Currently the following types of arrays are supported

"Sentrix96": Sentrix array, 12 columns, 8 rows

"Sentrix16": Sentrix array, 2 columns, 8 rows

"Slide12"  : Slide with 12 samples
}
\value{
  \code{arrayType} and \code{arrayID} return a \code{character} value
}
\author{Jan Oosting}
\keyword{ manip }
