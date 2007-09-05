\name{reportGenotypeSegmentation}
\alias{reportGenotypeSegmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reportGenotypeSegmentation(object, plotRaw = TRUE, subsample = NULL, panels = 0, minProbes = 10, maxY = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{plotRaw}{ ~~Describe \code{plotRaw} here~~ }
  \item{subsample}{ ~~Describe \code{subsample} here~~ }
  \item{panels}{ ~~Describe \code{panels} here~~ }
  \item{minProbes}{ ~~Describe \code{minProbes} here~~ }
  \item{maxY}{ ~~Describe \code{maxY} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function (object, plotRaw = TRUE, subsample = NULL, panels = 0, 
    minProbes = 10, maxY = 2, ...) 
{
    if (!all(c("lair", "nor.gt", "loh") \%in\% assayDataElementNames(object))) 
        stop("'calculateLOH' should be performed before making this report")
    if (!all(c("observed", "states", "predicted") \%in\% assayDataElementNames(object))) 
        stop("'segmentate' should be performed before making this report")
    subsample <- getSubsample(object, subsample)
    if (panels == 0) 
        panels <- length(levels(subsample))
    par(mfrow = c(panels, 1), mar = c(2.6, 4.1, 3.6, 1.6))
    for (smp in 1:ncol(object)) {
        for (subsmp in levels(subsample)) {
            selection <- subsample == subsmp
            chrs <- summary(as.factor(featureData(object)$CHR[selection]))
            selection <- selection & featureData(object)$CHR \%in\% 
                names(chrs[chrs >= 10])
            if (is.null(maxY) | maxY < 0) 
                maxY.opa <- max(assayData(object)$observed[selection, 
                  smp], na.rm = TRUE)
            else maxY.opa <- maxY
            plot(0, type = "n", main = sampleNames(object)[smp], 
                ylab = "intensity", yaxt = "none", xlab = "", 
                xaxt = "none", ...)
            par(usr = c(0, sum(selection), -0.25, maxY.opa))
            if (plotRaw) 
                points(assayData(object)$observed[selection, 
                  smp], pch = ".")
            chr <- featureData(object)$CHR[selection]
            xax <- getMidMaxIdx(chr)
            axis(1, xax$midpos, row.names(xax))
            axis(2)
            abline(v = xax$maxpos)
            points(assayData(object)$predicted[selection, smp], 
                pch = "-", col = "red")
            lair.offset <- -0.15
            lair.range <- 0.4
            abline(h = lair.offset + lair.range)
            het.nrm <- assayData(object)$nor.gt[selection, smp] == 
                "TRUE"
            het.nrm <- names(het.nrm)[het.nrm]
            het.nrm <- het.nrm[!is.na(het.nrm)]
            idx <- which(featureNames(object)[selection] \%in\% 
                het.nrm)
            points(idx, lair.offset + assayData(object)$lair[het.nrm, 
                smp] * lair.range, pch = "_", col = "blue")
            loh.offset <- -0.24
            loh.range <- 0.09
            loh.width <- 1.5
            q.col <- ifelse(assayData(object)$GSR[het.nrm, smp] < 
                0.8, "mediumblue", "green")
            col <- ifelse(assayData(object)$call[het.nrm, smp] == 
                "H", q.col, "red")
            segments(idx, loh.offset, idx, loh.offset + loh.range, 
                lwd = loh.width, col = col)
        }
    }
    invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
