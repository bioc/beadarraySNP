\name{copynumberConversion}
\alias{convert2aCGH}
\alias{convert2SegList}
\concept{SnpSetIllumina}

\title{Conversion to Copynumber analysis objects}

\description{
  SnpSetIllumina objects are converted to other objects for numerical analysis
}

\usage{
  convert2aCGH(object,normalizedTo=2,doLog=TRUE,organism="hsa")
  convert2SegList(object,normalizedTo=2,doLog=TRUE,organism="hsa")
}
\arguments{
  \item{object}{SnpSetIllumina object}
  \item{normalizedTo}{numeric, 'normal' copynumber datavalue for object}
  \item{doLog}{logical, perform logarithmic transformation (log2)}
  \item{organism}{character, organism used in object. Currently 'hsa' and
      'mmu' are recognized. Used to convert sex chromosomes to their 
      proper numerical representation}
}
\details{
  These functions produce objects that can be used by the analysis functions in
  the aCGH or snapCGH packages.
  The SnpSetIllumina intensity values are stored in a linear scale. Both types of
  objects assume a logarithmic scale, so by default the values are transformed to 
  a log2 scale centered around 0. 
}
\value{
  \code{convert2aCGH} returns a \code{aCGH} object as used in the aCGH package.
  \code{convert2SegList} returns a \code{SegList} object as used in the snapCGH package.
  
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}, \code{\link[aCGH:aCGH]{aCGH-class}}, \code{\link[snapCGH:SegList]{SegList-class}}
}
\keyword{manip}
