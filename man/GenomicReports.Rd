\name{GenomicReports}
\alias{reportChromosomesSmoothCopyNumber}
\alias{reportSamplesSmoothCopyNumber}
\alias{reportGenomeGainLossLOH}
\alias{reportChromosomeGainLossLOH}
\alias{reportGenomeIntensityPlot}
\title{Genomic reports}
\description{
  Create reports for all samples in a dataset.
}
\usage{
reportChromosomesSmoothCopyNumber(snpdata, grouping, normalizedTo=2, 
  smooth.lambda=2, ridge.kappa=0, plotLOH=c("none", "marker", "line", "NorTum"), 
  sample.colors = NULL, ideo.bleach=0.25, \dots)
reportSamplesSmoothCopyNumber(snpdata, grouping, normalizedTo=2, 
  smooth.lambda=2, ridge.kappa=0, plotLOH=c("none", "marker", "line", "NorTum"), 
  sample.colors=NULL, \dots)
reportGenomeGainLossLOH(snpdata, grouping, plotSampleNames=FALSE, sizeSampleNames=4, 
  distance.min, upcolor="red", downcolor="blue", lohcolor="grey", hetcolor="lightgrey", 
  lohwidth=1, segment=101, orientation=c("V","H"), \dots)
reportChromosomeGainLossLOH(snpdata, grouping, plotSampleNames=FALSE, distance.min,
  upcolor="red", downcolor="blue", lohcolor="grey", hetcolor="lightgrey", proportion=0.2, 
  plotLOH=TRUE, segment=101, \dots)
reportGenomeIntensityPlot(snpdata, normalizedTo=NULL, subsample=NULL, smoothing=c("mean", "quant"),
  dot.col="black", smooth.col="red", \dots)
}
\arguments{
  \item{snpdata}{SnpSetIllumina object.}
  \item{grouping}{factor, elements with same value are plotted together. Defaults
	  to groups of 4 in order of the samples in the object.}
  \item{normalizedTo}{numeric, a horizontal line is drawn at this position.}
  \item{smooth.lambda}{smoothing parameter for \code{\link[quantsmooth]{quantsmooth}}.}
  \item{ridge.kappa}{smoothing parameter for \code{\link[quantsmooth]{quantsmooth}}.}
  \item{plotLOH}{indicate regions or probes with LOH, see details.}
  \item{sample.colors}{vector of color.}
  \item{plotSampleNames}{logical.}
  \item{sizeSampleNames}{numeric, margin size for sample names.}
  \item{distance.min}{numerical.}
  \item{upcolor}{color.}
  \item{downcolor}{color.}
  \item{lohcolor}{color.}
  \item{hetcolor}{color.}
  \item{lohwidth}{}
  \item{segment}{integer.}
  \item{orientation}{["V","H"], vertical or horizontal orientation of plot.}
  \item{proportion}{}
  \item{subsample}{}
  \item{smoothing}{Type of smoothing per chromosome.}
  \item{dot.col}{color.}
  \item{smooth.col}{color.}
  \item{ideo.bleach}{numeric [0,1]}
  \item{\dots}{arguments are forwarded to \code{plot} or \code{getChangedRegions}.}
}
\details{
  The first function creates plots for each group and each chromosome in the 
  dataset. The second function creates full genome plot for each group in the
  dataset. Beware that a lot of plots can be created, and usually you should 
  prepare for that, by redirecting the plots to \code{pdf} or functions that
  create picture files like \code{jpg, png, bmp}.
}
\value{
  These functions are executed for their side effects
}
\author{Jan Oosting}
\seealso{
  \code{\link[quantsmooth]{quantsmooth}},\code{\link[quantsmooth]{prepareGenomeplot}},
  \code{\link{pdfChromosomesSmoothCopyNumber}}, \code{\link{pdfSamplesSmoothCopyNumber}}
}
\examples{
data(chr17.260)
chr17nrm <- standardNormalization(chr17.260)
par(mfrow = c(4,2), mar = c(2,4,2,1))
reportChromosomesSmoothCopyNumber(chr17nrm, grouping=pData(chr17.260)$Group,smooth.lambda = 4)
}
\keyword{manip}
