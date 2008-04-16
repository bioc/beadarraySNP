\name{smoothed.intensity}
\alias{smoothed.intensity}
\title{ Smooth intensity data }
\description{
  Create a table of smoothe intensity values
}
\usage{
smoothed.intensity(snpdata, smooth.lambda = 4, tau = 0.5)
}
\arguments{
  \item{snpdata}{SnpSetIllumina object}
  \item{smooth.lambda}{smoothing parameter}
  \item{tau}{quantile to smooth}
}
\value{
  Numerical matrix with same dimensions as data
}
\author{Jan Oosting}
\seealso{
   \code{\link{SnpSetIllumina-class}}
}
\keyword{manip}

