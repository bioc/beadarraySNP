\name{normalizeBetweenAlleles.SNP}
\alias{normalizeBetweenAlleles.SNP}
\concept{SnpSetIllumina}

\title{between Allele normalization}

\description{
   Perform between Allele normalization on Illumina Golden Gate bead arrays
}

\usage{
   normalizeBetweenAlleles.SNP(object,method=c("quantile"),subsample="OPA")
}
\arguments{
   \item{object}{\code{object} class SnpSetIllumina}
   \item{method}{\code{method} type of normalization}
   \item{subsample}{\code{subsample} factor or column name in \code{featureData} slot}
}

\details{
This function performs a quantile normalization between the Red and Green channels
for each sample. The rationale for this procedure stems from the fact that the
allele frequencies within each channel are always very similar, even in the 
presence of genomic abnormalities.
}
\value{
  This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}},\code{\link{normalizeWithinArrays.SNP}},
   \code{\link{backgroundCorrect.SNP}}
}

\examples{
data(chr17.260)
data.nrm<-normalizeBetweenAlleles.SNP(chr17.260)
}
\keyword{manip}
