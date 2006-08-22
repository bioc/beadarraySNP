\name{reportSamplePanelQC-methods}
\docType{methods}
\alias{reportSamplePanelQC-methods}
\alias{reportSamplePanelQC}
\alias{reportSamplePanelQC,QCIllumina-method}
\title{ reportSamplePanelQC  }
\description{
 Show raw intensity values for green and red channel for all samples in an 
 experiment
}  
\usage{
reportSamplePanelQC(object, by=10, legend=TRUE, ...)
}
\arguments{
  \item{object}{QCIllumina object}
  \item{by}{numeric, number of samples in each plot}
  \item{legend}{logical, create a final plot with a common legend for the
    barplots}
  \item{\dots}{arguments are forwarded to \code{barplot}}
}

\examples{
data(QC.260)
par(mfrow=c(2,2))
reportSamplePanelQC(QC.260,by=8)
}

\keyword{methods}
\keyword{hplot}
