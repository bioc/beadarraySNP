\name{GenomicReports}
\alias{reportChromosomesSmoothCopyNumber}
\alias{reportSamplesSmoothCopyNumber}
\title{Genomic reports}
\description{
  create reports for all samples in a dataset
}
\usage{
reportChromosomesSmoothCopyNumber(snpdata, grouping, normalizedTo = 2, 
  smooth.lambda = 2, ridge.kappa = 0, plotLOH = c("none", "marker", "line", "NorTum"), ...)
reportSamplesSmoothCopyNumber(snpdata, grouping, normalizedTo = 2, 
  smooth.lambda = 2, ridge.kappa = 0, plotLOH = c("none", "marker", "line", "NorTum"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{SnpSetIllumina object}
  \item{grouping}{factor, elements with same value are plotted together. Defaults
	  to groups of 4 in order of the samples in the object}
  \item{normalizedTo}{numeric, a horizontal line is drawn at this position}
  \item{smooth.lambda}{smoothing parameter for \code{\link[quantsmooth]{quantsmooth}}}
  \item{ridge.kappa}{smoothing parameter for \code{\link[quantsmooth]{quantsmooth}}}
  \item{plotLOH}{indicate regions or probes with LOH, see details}
  \item{\dots}{arguments are forwarded to \code{plot}}
}
\details{
  The first function creates plots for each group and each chromosome in the 
  dataset. The second function creates full genome plot for each group in the
  dataset. Beware that a lot of plots can be created, and usually you should 
  prepare for that, by redirecting the plots to \code{pdf} or functions that
  create picture files like \code{jpg, png, bmp}
}
\value{
}
\author{Jan Oosting}
\seealso{
  \code{\link[quantsmooth]{quantsmooth}},\code{\link[quantsmooth]{prepareGenomeplot}},
  \code{\link{pdfChromosomesSmoothCopyNumber}}, \code{\link{pdfSamplesSmoothCopyNumber}}
}
\examples{
data(chr17.260)
chr17nrm<-standardNormalization(chr17.260)
par(mfrow=c(4,2),mar=c(2,4,2,1))
reportChromosomesSmoothCopyNumber(chr17nrm, grouping=pData(chr17.260)$Group,smooth.lambda = 4)
}
\keyword{manip}
