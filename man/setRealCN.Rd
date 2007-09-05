\name{setRealCN}
\alias{setRealCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setRealCN(object, sample, cn.sum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{cn.sum}{ ~~Describe \code{cn.sum} here~~ }
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
function (object, sample, cn.sum) 
{
    if (!("inferred" \%in\% assayDataElementNames(object))) {
        inferred <- matrix(NA, nrow = nrow(object), ncol = ncol(object), 
            dimnames = dimnames(assayData(object)$G))
    }
    else {
        inferred <- assayData(object)$inferred
    }
    subsample <- beadarraySNP:::getSubsample(object, "OPA")
    for (subsmp in levels(subsample)) {
        selection <- subsample == subsmp
        st.sel <- which(cn.sum$states$opa == subsmp)
        st.pred <- cn.sum$states$intensity[st.sel]
        st.tcn <- cn.sum$states$copynumber[st.sel]
        trueCN <- st.tcn[match(assayData(object)$predicted[selection, 
            sample], st.pred)]
        inferred[selection, smp] <- trueCN
    }
    res <- object
    assayData(res)$inferred <- inferred
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
