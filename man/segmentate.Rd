\name{segmentate}
\alias{segmentate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segmentate(object, method = c("DNACopy", "HMM", "BioHMM", "GLAD"), normalizedTo = 2, doLog = TRUE, doMerge = FALSE, subsample = "OPA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{normalizedTo}{ ~~Describe \code{normalizedTo} here~~ }
  \item{doLog}{ ~~Describe \code{doLog} here~~ }
  \item{doMerge}{ ~~Describe \code{doMerge} here~~ }
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
function (object, method = c("DNACopy", "HMM", "BioHMM", "GLAD"), 
    normalizedTo = 2, doLog = TRUE, doMerge = FALSE, subsample = "OPA") 
{
    method <- match.arg(method)
    if (!require(snapCGH)) 
        stop("The snapCGH is not installed")
    sl <- convert2SegList(object, normalizedTo = normalizedTo, 
        doLog = doLog)
    subsample <- getSubsample(object, subsample)
    sl$genes <- cbind(sl$genes, subsample = subsample, pos.old = sl$genes[, 
        "Position"], chr.old = sl$genes[, "Chr"])
    sl$genes[, "Position"] <- sl$genes[, "Position"] + (sl$genes[, 
        "Chr"] * 1e+09)
    sl$genes[, "Chr"] <- rep(1, nrow(sl))
    subsamples <- unique(as.character(subsample))
    states <- matrix(NA, ncol = ncol(sl), nrow = nrow(sl), dimnames = dimnames(sl))
    predicted <- matrix(NA, ncol = ncol(sl), nrow = nrow(sl), 
        dimnames = dimnames(sl))
    for (smp in subsamples) {
        selection <- subsample == smp
        switch(method, DNACopy = {
            segm <- runDNAcopy(sl[selection, ])
        }, HMM = {
            segm <- runHomHMM(sorted.sl[selection, ], criteria = "AIC")
        }, BioHMM = {
            segm <- runBioHMM(sorted.sl[selection, ])
        }, GLAD = {
            segm <- runGLAD(sorted.sl[selection, ])
        })
        if (doMerge) 
            segm <- mergeStates(segm)
        states[selection, ] <- segm$state
        predicted[selection, ] <- segm$M.predicted
    }
    res <- object
    assayData(res)$observed <- sl$M.observed
    assayData(res)$states <- states
    assayData(res)$predicted <- predicted
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
