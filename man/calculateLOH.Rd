\name{calculateLOH}
\alias{calculateLOH}
\alias{calculateLair}
\title{Determine LOH in experiment}
\description{
  Using pairings of normal and tumor samples the LOH pattern is determined
}
\usage{
calculateLOH(object, grouping, NorTum = "NorTum", ...)
calculateLair(object, grouping = NULL, NorTum = "NorTum", min.intensity = NULL, 
    use.homozygous.avg = FALSE)
}
\arguments{
  \item{object}{SnpSetIllumina object}
  \item{grouping}{Factor to show which samples belong together (are of the 
    same individual)}
  \item{NorTum}{character vector or factor. Elements containing "N" are considered 
    to be the normal sample}
  \item{min.intensity}{numeric}
  \item{use.homozygous.avg}{logical}
  \item{\dots}{extra arguments for \code{link{heterozygousSNPs}}}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  For calculateLOH a SnpSetIllumina object with \code{loh} and \code{nor.gt} matrices in assayData. 
  \code{loh} is a logical matrix, and \code{nor.gt} is a character matrix containing
  the genotypes of the corresponding normal sample
  For calculateLair a SnpSetIllumina object with \code{lair} matrix in assayData. lair is the lesser
  allele intensity ratio. If a corresponding normal sample is found, it is taken as reference. Else the
  genotypes of normal samples are taken as a reference
}
\author{Jan Oosting}
\seealso{\code{\link{SnpSetIllumina-class}}}
\keyword{manip}
