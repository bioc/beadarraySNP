\name{plotGoldenGate4OPA}
\alias{plotGoldenGate4OPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotGoldenGate4OPA(object, cn.sum = NULL, sample = 1, plotRaw = FALSE, main = NULL, interact = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{cn.sum}{ ~~Describe \code{cn.sum} here~~ }
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{plotRaw}{ ~~Describe \code{plotRaw} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{interact}{ ~~Describe \code{interact} here~~ }
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
function (object, cn.sum = NULL, sample = 1, plotRaw = FALSE, 
    main = NULL, interact = FALSE, ...) 
{
    if (!all(c("lair", "nor.gt", "loh") \%in\% assayDataElementNames(object))) 
        stop("'calculateLOH' should be performed before making this plot")
    if (!all(c("observed", "states", "predicted") \%in\% assayDataElementNames(object))) 
        stop("'segmentate' should be performed before making this plot")
    subsample <- beadarraySNP:::getSubsample(object, "OPA")
    if (length(sample) != 1) 
        stop("Only 1 sample can be handled at a time")
    if (is.numeric(sample)) {
        if (sample < 1 | sample > ncol(object)) 
            stop(paste("sample", sample, "is not between 1 and ", 
                ncol(object)))
    }
    else {
        if (!(sample \%in\% sampleNames(object))) 
            stop(paste(sample, "is not a sample name in the object"))
    }
    if (is.null(cn.sum)) 
        cn.sum <- createCNSummary(object, 1)
    if (is.null(main)) {
        if (is.numeric(sample)) 
            main = sampleNames(object)[sample]
        else main = sample
    }
    par(mar = c(1, 3, 3, 1))
    plot(0, type = "n", main = main, ylab = "", yaxt = "none", 
        xlab = "", xaxt = "none")
    if (interact) 
        y.bottom <- 0
    else y.bottom <- 1.2
    par(usr = c(0, 1, y.bottom, 14))
    y.base <- 12
    for (subsmp in levels(subsample)) {
        selection <- subsample == subsmp
        chrs <- summary(as.factor(featureData(object)$CHR[selection]))
        selection <- selection & featureData(object)$CHR \%in\% 
            names(chrs[chrs >= 10])
        selected.n <- sum(selection)
        abline(h = y.base + c(-0.25, 0.25, 1, 2), col = c("black", 
            "black", "grey", "black"))
        if (plotRaw) 
            points((1:selected.n)/selected.n, y.base + assayData(object)$observed[selection, 
                sample], pch = ".")
        chr <- featureData(object)$CHR[selection]
        xax <- getMidMaxIdx(chr)
        text(xax$midpos/selected.n, y.base - 1, row.names(xax), 
            pos = 3)
        axis(2, y.base + c(0.5, 1, 1.5), c(0.5, 1, 1.5))
        segments(xax$maxpos/selected.n, y.base - 1, xax$maxpos/selected.n, 
            y.base + 2)
        predicted <- assayData(object)$predicted[selection, sample]
        points((1:selected.n)/selected.n, y.base + predicted, 
            pch = "-", col = "red")
        st.sel <- which(cn.sum$states$opa == subsmp)
        st.pred <- cn.sum$states$intensity[st.sel]
        st.tcn <- cn.sum$states$copynumber[st.sel]
        trueCN <- st.tcn[match(assayData(object)$predicted[selection, 
            sample], st.pred)]
        for (i in unique(trueCN)) {
            probes <- trueCN == i
            predicted[probes] <- mean(predicted[probes])
        }
        points((1:selected.n)/selected.n, y.base + predicted, 
            pch = "-", col = "green")
        idx <- which(c(trueCN, trueCN[length(trueCN)]) != c(-100, 
            trueCN) | c(chr[1], chr) != c(chr, chr[length(chr)]))
        posi <- (idx + c(idx[-1], length(trueCN) + 1))/2 - 0.5
        text(posi/selected.n, y.base - 0.7, trueCN[idx], pos = 3, 
            col = "green")
        segments(idx/selected.n, y.base - 1, idx/selected.n, 
            y.base - 0.25, col = "green")
        lair.offset <- -0.15
        lair.range <- 0.4
        het.nrm <- assayData(object)$nor.gt[selection, sample] == 
            "TRUE"
        het.nrm <- names(het.nrm)[het.nrm]
        het.nrm <- het.nrm[!is.na(het.nrm)]
        idx <- which(featureNames(object)[selection] \%in\% het.nrm)
        points(idx/selected.n, y.base + lair.offset + assayData(object)$lair[het.nrm, 
            sample] * lair.range, pch = "-", col = "blue")
        loh.offset <- -0.25
        loh.range <- 0.1
        loh.width <- 1.5
        q.col <- ifelse(assayData(object)$GSR[het.nrm, sample] < 
            0.8, "mediumblue", "green")
        col <- ifelse(assayData(object)$call[het.nrm, sample] == 
            "H", q.col, "red")
        segments(idx/selected.n, y.base + loh.offset, idx/selected.n, 
            y.base + loh.offset + loh.range, lwd = loh.width, 
            col = col)
        y.base <- y.base - 3
    }
    abline(h = 2)
    if (interact) {
        rect(c(0.9, 0.9), c(0, 1), c(1, 1), c(1, 2))
        text(0.9, c(0.5, 1.5), c("Redraw", "Raw data"), pos = 4)
    }
    text(0.02, y.bottom + 0.4, paste("Target Index", cn.sum$dnaIndex, 
        "Current index", format(getDNAindex(cn.sum), digits = 2, 
            nsmall = 2)), adj = 0)
    invisible(cn.sum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
