\name{backgroundCorrect.SNP}
\alias{backgroundCorrect.SNP}
\concept{SnpSetIllumina}

\title{Background correction}

\description{
   Perform background correction on Illumina Golden Gate bead arrays
}

\usage{
   backgroundCorrect.SNP(object,method=c("none", "subtract", "half", "minimum",
                         "edwards", "normexp", "rma"),offset = 0)
}
\arguments{
   \item{object}{SnpSetIllumina object}
   \item{method}{character, method of correction}
   \item{offset}{numeric, constant to add after correction}
}

\details{
  Code has been ported from the \code{limma} package. The matrices \code{Gb} and
  \code{Rb} should be available in the \code{arrayData} slot of the object.
}
\value{
  This function returns an \code{SnpSetIllumina} object with background 
  corrected values in the \code{G} and \code{R}.
}
\author{Jan Oosting, based on \code{limma} package by G. Smyth}

\seealso{
   \code{\link{SnpSetIllumina-class}}, \code{\link[limma]{backgroundCorrect}},
   \code{\link{backgroundEstimate}}, \code{\link{normalizeBetweenAlleles.SNP}},
   \code{\link{normalizeWithinArrays.SNP}}
}

\examples{
  \dontrun{data.bg<-backgroundCorrect.SNP(data.raw,"subtract")}
}
\keyword{manip}
