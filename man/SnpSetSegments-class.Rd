\name{SnpSetSegments-class}
\Rdversion{1.1}
\docType{class}
\alias{SnpSetSegments-class}
\alias{cn.segments,SnpSetSegments-method}
\alias{cn.segments<-,SnpSetSegments,list-method}
\alias{cn.segments}
\alias{cn.segments<-}
\alias{initialize,SnpSetSegments-method}

\title{Class "SnpSetSegments"}
\description{The SnpSetSegments class is a direct descendant of the SnpSetIllumina class, with an extra slot to define the genomic segments in each sample. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpSetSegments", assayData, phenoData, experimentData, annotation, protocolData, call, callProbability, G, R, cn.segments, featureData, extraData, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{cn.segments}:}{Object of class \code{"list"}}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{annotation}:}{Object of class \code{"character"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} see \code{"\linkS4class{SnpSetIllumina}"} }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} \code{"\linkS4class{VersionedBiobase}"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpSetIllumina}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "SnpSetIllumina", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpSetIllumina", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "SnpSetIllumina", distance 4.
}
\section{Methods}{
  \describe{
    \item{cn.segments}{\code{signature(object = "SnpSetSegments")}: ... }
    \item{cn.segments<-}{\code{signature(object = "SnpSetSegments", value = "list")}: ... }
    \item{initialize}{\code{signature(.Object = "SnpSetSegments")}: ... }
	 }
}
\references{ Corver et.al. Can Res dec 2008 }
\author{ Jan Oosting }
\note{ This class is under development, and not usable in the current form }

\seealso{
	\code{\link{segmentate}}
}
\examples{
showClass("SnpSetSegments")
}
\keyword{classes}
