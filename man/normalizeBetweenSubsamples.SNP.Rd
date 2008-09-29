\name{normalizeBetweenSubsamples.SNP}
\alias{normalizeBetweenSubsamples.SNP}
\title{ Normalization between subsamples }
\description{
  Quantile normalization between subsamples within a single SnpSetIllumina object
}
\usage{
normalizeBetweenSubsamples.SNP(object, subsample = "OPA")
}
\arguments{
  \item{object}{class \code{SnpSetIllumina}}
  \item{subsample}{\code{factor} with length number of features in object or \code{char},column name in \code{featureData} slot}
}
\details{
  Perform quantile normalization of the red and green channel between subsamples. This can be used in situations where multiple different 
  assays that cover the same genomic regions (or whole genome) have been done on the same biological specimen. 
  This function was introduced for version 5 Golden Gate Linkage analysis that consist of 4 assays of ~ 1500 probes. Where previous versions 
  of this assay each targeted a number of chromosomes, in version 5 each assay covers the whole genome.
}
\value{
  This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}},\code{\link{normalizeBetweenAlleles.SNP}},
   \code{\link{normalizeWithinArrays.SNP}},\code{\link{backgroundCorrect.SNP}}
}
\examples{
data(chr17.260)
data.nrm<-normalizeBetweenSubsamples.SNP(chr17.260)
}
\keyword{manip}
