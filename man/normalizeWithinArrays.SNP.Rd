\name{normalizeWithinArrays.SNP}
\alias{normalizeWithinArrays.SNP}
\concept{SnpSetIllumina}

\title{within Array normalization}

\description{
   Perform within array normalization on Illumina Golden Gate bead arrays
}

\usage{
   normalizeWithinArrays.SNP(object, callscore=0.5, normprob=0.5, quantilepersample=FALSE,
                             relative=FALSE, fixed=FALSE, useAll=FALSE, subsample="OPA")
}
\arguments{
   \item{object}{class SnpSetIllumina}
   \item{callscore}{numeric with range 0:1, threshold for probe inclusion}
   \item{normprob}{numeric with range 0:1, target quantile for normalization. 
     The default is to divide the sample intensities by the median of the 
     selected subset}
   \item{quantilepersample}{logical. If \code{TRUE} then the threshold is 
     determined for each sample, else it is experiment wide. This is only 
     relevant when \code{fixed} is \code{FALSE}. }
   \item{relative}{logical. If \code{TRUE} then the ratio of GCS and GTS is used,
     else only the GCS is used as the quality score}
   \item{fixed}{logical. If \code{TRUE} then \code{callscore} is the fixed 
     threshold for the quality score, else the probes above the quantile 
     \code{callscore} are used.}
   \item{useAll}{logical. If \code{TRUE} then all probes in the dataset are 
     eligible as the invariant set, else only the heterozygous SNPs}
   \item{subsample}{factor or column name in \code{featureData} slot, the levels 
    of the factor are treated separately}
}

\details{
  The function uses high quality heterozygous SNPs as an invariant set with the
assumption that these have the highest probability of coming from 
unaffected regions of the genome. Most of the arguments are used to determine
the quality of the call.\cr
 
}
\value{
  This function returns a \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina}},\code{\link{normalizeLoci.SNP}},
   \code{\link{backgroundCorrect.SNP}},\code{\link{normalizeBetweenAlleles.SNP}}
}

\examples{
  data(chr17.260)
  data.nrm<-normalizeWithinArrays.SNP(chr17.260)

}
\keyword{manip}
