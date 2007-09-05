\name{getIdGoldenGate4OPA}
\alias{getIdGoldenGate4OPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getIdGoldenGate4OPA(object, opa, xpos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{opa}{ ~~Describe \code{opa} here~~ }
  \item{xpos}{ ~~Describe \code{xpos} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, opa, xpos) 
{
    subsample <- beadarraySNP:::getSubsample(object, "OPA")
    subsmp <- levels(subsample)[opa]
    selection <- subsample == subsmp
    chrs <- summary(as.factor(featureData(object)$CHR[selection]))
    selection <- featureNames(object)[selection & featureData(object)$CHR \%in\% 
        names(chrs[chrs >= 10])]
    selection[ceiling(xpos * length(selection))]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
