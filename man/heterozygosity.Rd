\name{heterozygosity}
\alias{heterozygosity}
\title{Find regions of homozygous SNPs}
\description{
  Analyze affected material without corresponding unaffected material in order
  to find regions that contain stretches of homozygous SNPs as an indication of
  loss of heterozygosity (LOH)
}
\usage{
heterozygosity(genotype, decay = 0.8, threshold = 0.1)
}
\arguments{
  \item{genotype}{character or logical vector, genotypes of affected tissue}
  \item{decay}{numeric in range (0,1)}
  \item{threshold}{numeric in range (0,1)}
}
\details{
  The method calculates how long the stretch of homozygous SNPs is for each element
  \code{delay} and \code{threshold} can be set to skip individual heterozygous 
  probes in a longer stretch of homozygous probes. The default setting 
  tolerate 1 erroneous heterozygous SNP every 10 homozygous SNPs. Set \code{threshold} 
  at 1 to stop discarding hetrozygous SNPs
}
\value{
	A numeric vector with the same length as genotype is returned. Higher values,
	of 15 and higher, indicate regions of LOH
}
\author{Jan Oosting}
\seealso{ \code{\link{compareGenotypes}}, \code{\link{heterozygousSNPs}}}
\examples{
  data(chr17.260)
  plot(heterozygosity(exprs(chr17.260)[,"514TV"]))
}
\keyword{manip}
