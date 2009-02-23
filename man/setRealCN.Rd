\name{setRealCN}
\alias{setRealCN}
\title{Integrate state information into SNP object}
\description{
  Set calculated values of copy numbers in \code{inferred} element of \code{AssayData} slot
}
\usage{
setRealCN(object, sample, cn.sum, subsample="OPA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{class SnpSetIllumina }
  \item{sample}{sample identifier}
  \item{cn.sum}{list, see \code{createCNSummary}}
  \item{subsample}{"OPA"}
}
\value{
  SnpSetIllumina object with \code{inferred} element of \code{AssayData} slot set 
}

\author{Jan Oosting}
\seealso{ \code{\link{segmentate}, \link{alterCN}, \link{plotGoldenGate4OPA} \link{createCNSummary}} }
\keyword{manip}
