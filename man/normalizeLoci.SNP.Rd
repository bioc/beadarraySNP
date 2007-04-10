\name{normalizeLoci.SNP}
\alias{normalizeLoci.SNP}
\concept{SnpSetIllumina}

\title{locus normalization}

\description{
   Perform locus normalization on Illumina Golden Gate bead arrays
}

\usage{
   normalizeLoci.SNP(object, method=c("normals","paired"), NorTum="NorTum",
                     Gender="Gender", Subject="Subject", normalizeTo=2, trig=FALSE)
}
\arguments{
   \item{object}{\code{object} class SnpSetIllumina}
   \item{method}{character. If \code{"normals"} then all normal samples in the 
     dataset are used as the invariant set. If \code{"paired"} then affected 
     samples are normalized to their paired normal samples}
   \item{NorTum}{logical or character vector or name of column in \code{pData} slot. 
     depicts the normal, unaffected samples in the dataset. In a character vector
     these should have the value \code{"N"}}
   \item{Gender}{logical or character vector or name of column in \code{pData} slot.
     depicts the female samples in the dataset and is used to normalize the sex
     chromosomes. In a character vector these should have value \code{"F"}}
   \item{Subject}{factor or name of or column in \code{pData} slot. This factor
     is used to pair the samples when \code{method} is \code{"paired"}}
   \item{normalizeTo}{\code{normalizeTo} numeric. The average copy number of the
     sample.}
   \item{trig}{Logical, use geometric distance of intensity. Otherwise use 
     addition of intensities}
}

\details{
  This function is usually performed in the last step of normalization in order
  to obtain calculated copy numbers.
  
}
\value{
  This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina}},\code{\link{normalizeWithinArrays.SNP}},
   \code{\link{normalizeBetweenAlleles.SNP}}
}

\examples{
  data(chr17.260)
  data.nrm<-normalizeLoci.SNP(chr17.260)

}
\keyword{manip}
