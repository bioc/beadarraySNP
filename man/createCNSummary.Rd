\name{createCNSummary}
\alias{createCNSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createCNSummary(object, sample, dnaIndex, subsample = "OPA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{dnaIndex}{ ~~Describe \code{dnaIndex} here~~ }
  \item{subsample}{ ~~Describe \code{subsample} here~~ }
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
function (object, sample, dnaIndex, subsample = "OPA") 
{
    deftarget <- round(dnaIndex * 2)
    subsample <- beadarraySNP:::getSubsample(object, subsample)
    res <- NULL
    getExisting <- FALSE
    if (!is.null(assayData(object)$inferred)) 
        if (!any(is.na(assayData(object)$inferred[, sample]))) 
            getExisting <- TRUE
    for (subsmp in levels(subsample)) {
        selection <- subsample == subsmp
        tmp <- summary(as.factor(assayData(object)$predicted[selection, 
            sample]))
        if (getExisting) 
            deftarget <- aggregate(assayData(object)$inferred[selection, 
                sample], by = list(assayData(object)$predicted[selection, 
                sample]), FUN = mean)[, 2]
        res <- rbind(res, data.frame(opa = rep(subsmp, length(tmp)), 
            count = tmp, intensity = names(tmp), copynumber = deftarget))
    }
    copynumber.total.nrm <- sum(numericCHR(featureData(object)$CHR) < 
        90) * 2
    if (!is.null(pData(object)$Gender)) 
        gender <- pData(object)[sample, "Gender"]
    else gender <- "F"
    if (gender == "F") 
        copynumber.total.nrm <- copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR) == 
            98) * 2
    else copynumber.total.nrm <- copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR) \%in\% 
        98:99)
    list(dnaIndex = dnaIndex, CN.total.nrm = copynumber.total.nrm, 
        states = res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
