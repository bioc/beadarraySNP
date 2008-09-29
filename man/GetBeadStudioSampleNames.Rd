\name{GetBeadStudioSampleNames}
\alias{GetBeadStudioSampleNames}
\title{Extract samplenames from a report file}
\description{
  Extract the samplenames from a report file that was created as a final report from Illumina Beadstudio
}
\usage{
GetBeadStudioSampleNames(reportfile)
}
\arguments{
  \item{reportfile}{character, name of report file }
}
\details{
This function will read the report file, and extract the sample names from the \code{Sample ID} column
}
\value{
  character vector
}
\author{Jan Oosting}
\seealso{ \code{\link{read.SnpSetIllumina}}}
\keyword{manip}
