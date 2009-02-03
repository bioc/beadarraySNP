\name{backgroundEstimate}
\alias{backgroundEstimate}
\concept{SnpSetIllumina}

\title{Estimate background intensities from foreground intensity}

\description{
  Background intensity from Illumina Golden Gate bead arrays are estimated
  based on several data models
}

\usage{
  backgroundEstimate(object,method=c("minimum", "mode","intmin", 
  "anglemode"), maxmode=3000, bincount=40, maxangle=0.3, subsample="OPA")
}
\arguments{
  \item{object}{SnpSetIllumina object}
  \item{method}{chracter, data model to use}
  \item{maxmode}{numeric, maximum intensity for mode for \code{method="mode"}}
  \item{bincount}{numeric, for \code{method="intmin"} , see details}
  \item{maxangle}{numeric in radians, maximum theta for mode for \code{method="anglemode"} }
  \item{subsample}{factor or column name in \code{featureData} slot}
}

\details{
  The Illumina software does not provide background values in the output. 
  Some models can be used to estimate background from the raw data intensities.\cr
  \code{minimum}: The allele specific minimum intensity is used.\cr
  \code{mode}: This model assumes that the first mode of the density of the 
    intensities is determined by the zero-allele in the data, see ref. The 
    signal intensity of the zero-allele should be zero, therefore this is 
    considered the background value.\cr
  \code{intmin}: This model assumes there is crosstalk between the alleles, and 
    background increases with the intensity of the other allele. The range 
    between 0 and the maximum of the other allele is divided in \code{bincount}
    bins, and the minimum for this allele is determined for probes where the 
    other allele falls in a bin. A linear fit is determined though the minimum
    values to obtain a gradually increasing value.\cr
  \code{anglemode}: This model finds the density modes closest to 0 and \eqn{\frac{\pi}{2}}  
    for polar transformed intensities, and uses this to determine background.
}
\value{
  This function returns an \code{SnpSetIllumina} object. The \code{Rb} and \code{Gb} 
  matrices in the \code{assayData} slot contain estimated background values.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}, \code{\link{backgroundCorrect.SNP}}
}
\keyword{manip}
