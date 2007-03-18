\name{copynumberConversion}
\alias{convert2aCGH}
\alias{convert2SegList}
\concept{SnpSetIllumina}

\title{Conversion to Copynumber analysis objects}

\description{
  SnpSetIllumina objects are converted to other objects for numerical analysis
}

\usage{
  convert2aCGH(object,normalizedTo=2,organism="hsa")
  convert2SegList(object,normalizedTo=2,organism="hsa")
}
\arguments{
  \item{object}{SnpSetIllumina object}
  \item{normalizedTo}{numeric, 'normal' copynumber datavalue for object}
  \item{organism}{character, organism used in object. Currently 'hsa' and
      'mmu' are recognized. Used to convert sex chromosomes to their 
      proper numerical representation}
}
\details{
  These functions produce objects that can be used by the analysis functions in
  the aCGH or snapCGH packages.
  The SnpSetIllumina intensity values are stored in a linear scale. Both types of
  objects assume a logarithmic scale, so the values are transformed to a log2 scale 
  centered around 0. 
}
\value{
  \code{convert2aCGH} returns a \code{aCGH} object as used in the aCGH package.
  \code{convert2SegList} returns a \code{SegList} object as used in the snapCGH package.
  This function returns an \code{SnpSetIllumina} object. The \code{Rb} and \code{Gb} 
  matrices in the \code{assayData} slot contain estimated background values.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}
}
\examples{

}
\keyword{manip}
