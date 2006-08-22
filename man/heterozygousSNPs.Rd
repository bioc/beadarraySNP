\name{heterozygousSNPs}
\alias{heterozygousSNPs}
\concept{SnpSetIllumina}

\title{Retrieve heterozygous SNPs}

\description{
   Heterozygous SNPs are determined based on quality score criteria
}

\usage{
  heterozygousSNPs(object, threshold=0.9, useQuality=TRUE, relative=TRUE, 
                   percentile=FALSE)}
\arguments{
   \item{object}{class SnpSetIllumina}
   \item{threshold}{numeric (0:1) minimum quality score to be called heterozygous}
   \item{useQuality}{logical, use quality score}
   \item{relative}{logical, use quality score relative to GTS, see details}
   \item{percentile}{logical, use percentage of probes above threshold}
}

\details{
  This function presumes that the specificity for determining heterozygity is 
more important than the sensitivity, and will therefore only call probes heterozygous 
if that can be done with high certainty.\cr
  The Illumina genotyping software calculates two quality measures: gen train score (GTS) 
and gen call score (GCS). The GTS is a measure for how well clusters can be recognized 
in a training set. This value is probe specific, and the same for all samples in an
experiment. The GCS is a probe-specific, sample specific value that measures how close 
a probe in a sample is to the clusters determined in the training step. This value is 
always lower than the GTS for a probe.\cr
  \code{read.SnpSetIllumina} will put GCS into the \code{callProbability} element 
of the \code{assaydata} slot and the GTS into the \code{reporterInfo} slot. The 
function uses these locations to retrieve the necessary information.\cr
  If \code{relative} is \code{FALSE} then the raw GCS values are compared to the
\code{threshold}. In this case a \code{threshold} of around 0.5 should be used. If 
\code{relative} is \code{TRUE} then GCS/GTS is compared to the \code{threshold} and
\code{threshold} should be around 0.9.\cr
  With \code{percentile=TRUE} the \code{threshold} quantile is calculated for each sample,
and only probes with higher scores can be called heterozygous. A \code{threshold}
of around 0.2 seems to work fine usually.
}
\note{
The purpose of the function is to separate heterozygous probes from non-heterozygous
probes. In tumor samples the determination of the genotype can be difficult, 
because of aneuploidy and the fact that a sample is often a mixture
of normal and tumor cells. 
}

\value{
  This function returns a \code{logical matrix} with same dimensions as \code{object}.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}}
}
\examples{
  data(chr17.260)
  plot(heterozygosity(heterozygousSNPs(chr17.260[,"514TV"])),col="red",pch="x")
  points(heterozygosity(exprs(chr17.260)[,"514TV"]))
}
\keyword{manip}
