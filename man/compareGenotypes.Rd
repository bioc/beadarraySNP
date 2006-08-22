\name{compareGenotypes}
\alias{compareGenotypes}
\title{ Compare genotypes }
\description{
  Pairwise comparison of genotypes between unaffected and affected tissue from the same subject
}
\usage{
compareGenotypes(genotypeT, genotypeN)
}
\arguments{
  \item{genotypeT}{character or logical vector, genotypes of affected tissue}
  \item{genotypeN}{character or logical vector with same length as genotypeT, genotypes of 
	  unaffected, normal tissue}
}
\details{
  Heterozygous probes have one the following values. \code{TRUE}, \code{'H'} or \code{'AB'}. 
	All other values are considered homozygous. The primary purpose of the method 
	is to find probes with loss of heterozygosity (LOH), where the unaffected probe 
	is heterozygous and the affected is called homozygous.
}
\value{
	A vector with the same length as the arguments where each element can have one of four values
	\item{'u'}{Uninformative: both affected and normal are homozygous}
	\item{'i'}{Informative: both affected and unaffected heterozygous}
  \item{'l'}{Loss: unaffected heterozygous, affected homozygous}
  \item{'a'}{Artefact: unaffected homozygous, affected heterozygous}
}
\author{Jan Oosting}
\seealso{ \code{\link{heterozygousSNPs}} }
\examples{
data(chr17.260)
compareGenotypes(exprs(chr17.260)[,"514TV"],exprs(chr17.260)[,"514NP"])
}
\keyword{manip}
