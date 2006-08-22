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
\value{
  \code{arrayType} and \code{arrayID} return a \code{character} value
}
\author{Jan Oosting}
\keyword{ manip }
