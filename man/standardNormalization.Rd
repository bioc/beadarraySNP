\name{standardNormalization}
\alias{standardNormalization}
\title{Default complete normalization}
\description{
  Performs all steps in normalization at best settings as determined in ref.
}
\usage{
standardNormalization(snpdata)
}
\arguments{
  \item{snpdata}{SnpSetIllumina object with raw data}
}
\details{
 The function performs in the following steps
 \code{snpdata<-normalizeBetweenAlleles.SNP(snpdata)}\cr
 \code{snpdata<-normalizeWithinArrays.SNP(snpdata,callscore = 0.8, relative = TRUE, 
       fixed = FALSE, quantilepersample = TRUE)}\cr
 \code{snpdata<-normalizeLoci.SNP(snpdata,normalizeTo = 2)}
 
}
\value{
  A \code{SnpSetIllumina} object with the \code{G, R and intensity} elements in 
  assayData normalized to obtain values close to 2 on a linear scale for 
  unaffected material.
}
\references{}
\author{Jan Oosting}
\seealso{\code{\link{normalizeBetweenAlleles.SNP}},\code{\link{normalizeWithinArrays.SNP}},
  \code{\link{normalizeLoci.SNP}}}
\examples{
  data(chr17.260)
  data.nrm<-standardNormalization(chr17.260)
}
\keyword{manip}
