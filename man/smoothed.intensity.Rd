\name{smoothed.intensity}
\alias{smoothed.intensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smooth intensity data }
\description{
  Create a table of smoothe intensity values
}
\usage{
smoothed.intensity(snpdata, smooth.lambda = 4, tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
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
   \code{\link{SnpSetIllumina-class}}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

