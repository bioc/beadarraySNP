\name{dist.GT}
\alias{dist.GT}
\title{ dist.GT }
\description{
  Calculate distance matrix based of differences in genotype calls
}
\usage{
dist.GT(object)
}
\arguments{
  \item{object}{SnpSetIllumina object}
}
\details{
  Calculates distances between samples as percentage of differences in genotype
}
\value{
  'dist.GT' returns an object of class 'dist'
}
\author{Jan Oosting}
\seealso{ \code{\link{dist}}, \code{\link{hclust}} }
\examples{
  data(chr17.260)
  plot(hclust(dist.GT(chr17.260)))
}
\keyword{ manip }
