\name{createCNSummary}
\alias{createCNSummary}

\title{ Summarization of Copy number states }
\description{
  Create a summary object of the genomic copy number states in a sample of segmented data
}
\usage{
createCNSummary(object, sample, dnaIndex=1, subsample = "OPA")
}

\arguments{
  \item{object}{ SNPSetIllumina object after segmentation \code{\link{segmentate}} }
  \item{sample}{ SampleName or index of the sample for which to create the summary}
  \item{dnaIndex}{ Measured DNA index of the sample}
  \item{subsample}{ factor or column name in \code{featureData} slot }
}
\details{
  The segments within a sample are assigned a copy number value. When the \code{inferred} slot 
  in \code{assayData} is empty, all segments will be set to 2. Otherwise the values are recovered
  from the \code{inferred} slot. Gender is taken into account for the sex chromosomes.  
}
\value{
  list with the following elements
  \item{dnaIndex }{same as parameter dnaIndex}
  \item{CN.total.nrm}{Total expected copynumber for a 'normal' specimen ~ 2*featurecount}
  \item{states}{data.frame with columns \code{opa, count ,intensity, copynumber}}
  This list can be used as the \code{cn.sum} argument for \code{plotGoldenGate4OPA, alterCN, 
  getDNAindex} and \code{setRealCN}
}
\author{ Jan Oosting}
\seealso{ \code{\link{segmentate}, \link{alterCN}, \link{plotGoldenGate4OPA} } }
\keyword{manip}
