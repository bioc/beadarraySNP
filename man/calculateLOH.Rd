\name{calculateLOH}
\alias{calculateLOH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine LOH in experiment}
\description{
  Using pairings of normal and tumor samples the LOH pattern is determined
}
\usage{
calculateLOH(object, grouping, NorTum = "NorTum", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{SnpSetIllumina object}
  \item{grouping}{Factor to show which samples belong together (are of the 
    same individual)}
  \item{NorTum}{character vector or factor. Elements containing "N" are considered 
    to be the normal sample}
  \item{\dots}{extra arguments for \code{link{heterozygousSNPs}}}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  SnpSetIllumina object with \code{loh} and \code{nor.gt} matrices in assayData. 
  \code{loh} is a logical matrix, and \code{nor.gt} is a character matrix containing
  the genotypes of the corresponding normal sample
}
\author{Jan Oosting}
\seealso{\code{\link{SnpSetIllumina-class}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
