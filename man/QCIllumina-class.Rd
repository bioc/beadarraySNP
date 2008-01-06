\name{QCIllumina-class}
\docType{class}
\alias{QCIllumina}
\alias{QCIllumina-class}
\alias{arrayID,QCIllumina-method}
\alias{arrayID<-,QCIllumina-method}
\alias{arrayType,QCIllumina-method}
\alias{arrayType<-,QCIllumina-method}
\alias{initialize,QCIllumina-method}
\alias{plotQC,QCIllumina-method}

\title{Class "QCIllumina" }
\description{
  Container of QC information on arrays that contain multiple samples. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QCIllumina", arrayType, arrayID, 
intensityMed, greenMed, redMed, intensityMode, greenMode, redMode, validn, annotation, 
samples)}, but are usually created by \code{\link{calculateQCarray}}.
}

\section{Slots}{
	 \describe{
    \item{\code{arrayType}:}{character, Type of array. See  \code{\link{arrayType}} }
    \item{\code{arrayID}:}{character, Array ID }
    \item{\code{intensityMed}:}{numeric matrix, Median of intensity of samples }
    \item{\code{greenMed}:}{numeric matrix, Median of green values }
    \item{\code{redMed}:}{numeric matrix, Median of red values }
    \item{\code{validn}:}{numeric matrix, Number of valid probe values in samples }
    \item{\code{annotation}:}{character matrix, Annotation of samples }
    \item{\code{samples}:}{character matrix, Sample IDs }
  }
}
\section{Methods}{
  \describe{
    \item{arrayID}{\code{signature(object = "QCIllumina")}: Returns type of array }
    \item{arrayID<-}{\code{signature(object = "QCIllumina")}: Sets type of array. 
		      Currently only \code{"Sentrix"} is supported }
    \item{arrayType}{\code{signature(object = "QCIllumina")}: Returns ID of array }
    \item{arrayType<-}{\code{signature(object = "QCIllumina")}: Sets ID/Barcode of array }
    \item{initialize}{\code{signature(.Object = "QCIllumina")}  }
    \item{plotQC}{\code{signature(object = "QCIllumina")}\code{character}: plots spatial
		     overview of QC information, \code{type} is one of \code{c("intensityMed",
				 "greenMed","redMed","validn","annotation","samples")}   }
	 }
}
\author{Jan Oosting}
\seealso{ \code{\link{calculateQCarray}} }
\keyword{classes}
