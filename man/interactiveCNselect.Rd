\name{interactiveCNselect}
\alias{interactiveCNselect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interactiveCNselect(object, sample = 1, dnaIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{dnaIndex}{ ~~Describe \code{dnaIndex} here~~ }
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
function (object, sample = 1, dnaIndex) 
{
    plotRaw <- TRUE
    cn.res <- createCNSummary(object, sample = sample, dnaIndex = dnaIndex)
    cn.res <- plotGoldenGate4OPA(object, cn.res, sample = sample, 
        plotRaw = plotRaw, interact = TRUE)
    repeat {
        location <- locator(n = 1)
        if (is.null(location)) 
            break
        if (location$x > 0.9 & location$y > 0 & location$y < 
            1) {
            cn.res <- plotGoldenGate4OPA(object, cn.res, sample = sample, 
                plotRaw = plotRaw, interact = TRUE)
        }
        else if (location$x > 0.9 & location$y > 1 & location$y < 
            2) {
            plotRaw <- !plotRaw
            cn.res <- plotGoldenGate4OPA(object, cn.res, sample = sample, 
                plotRaw = plotRaw, interact = TRUE)
        }
        else {
            opa <- floor((17 - location$y)/3)
            if (opa \%in\% 1:4) {
                snpid <- getIdGoldenGate4OPA(object, opa, location$x)
                increase <- ((16 - location$y) - opa * 3) < 0
                cn.res <- alterCN(cn.res, opa, assayData(object)$predicted[snpid, 
                  sample], increase)
                cn.res <- plotGoldenGate4OPA(object, cn.res, 
                  sample = sample, plotRaw = plotRaw, interact = TRUE)
                text(0, 1.6, paste("Chr", reporterInfo(object)[snpid, 
                  "CHR"], "position", reporterInfo(object)[snpid, 
                  "MapInfo"]), adj = 0)
                text(0, 1.2, paste("CN", assayData(object)$predicted[snpid, 
                  sample]), adj = 0)
                text(0, 0.8, paste("increase", increase), adj = 0)
            }
        }
    }
    cn.res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
