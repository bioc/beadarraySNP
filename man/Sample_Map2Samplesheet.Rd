\name{Sample_Map2Samplesheet}
\alias{Sample_Map2Samplesheet}
\title{Convert Beadstudio Sample Map file to samplesheet}
\description{
  Create a samplesheet that can be used to import Illumina beadstudio data
}
\usage{
Sample_Map2Samplesheet(samplemapfile, saveas = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplemapfile}{character, name of the SampleMap file}
  \item{saveas}{character, optional, name of samplesheet file that can be used
    directly by \code{\link{read.SnpSetIllumina}}}
}
\details{
  During the creation of a final reportfile from Beadstudio there is an option to create
  Map files. The \code{Sample_Map.txt} file can be used to create an initial samplesheet
  for use in the \code{read.SnpSetIllumina} function
}
\value{
  A \code{data.frame} with the samplesheet
}
\author{J. Oosting}
\seealso{\code{\link{read.SnpSetIllumina}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
