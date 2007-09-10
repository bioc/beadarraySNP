\name{getDNAindex}
\alias{getDNAindex}
\title{ Calculate the DNA index based on assigned copy number values to probes}
\description{
  Calculate the DNA index based on assigned copy number values to probes
}
\usage{
getDNAindex(cn.sum)
}
\arguments{
  \item{cn.sum}{list with elements \code{dnaIndex, CN.total.nrm, states}, see \code{\link{createCNSummary}}}
}
\details{
}
\value{
  scalar. DNA index of an unaffected sample is 1
}
\author{Jan Oosting}
\seealso{\code{\link{createCNSummary}, \link{plotGoldenGate4OPA}}}
\keyword{manip}
