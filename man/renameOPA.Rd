\name{renameOPA}
\alias{renameOPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Change the linkage panel in a dataset}
\description{
  Change the linkage panel in a dataset
}
\usage{
renameOPA(snpdata, newOPA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{SnpSetIllumina object}
  \item{newOPA}{character, new linkage panel}
}
\details{
  In order to combine different versions of the linkage panels, this function makes it 
  possible to map the equal SNPs in both datasets.
}
\value{
  SnpSetIllumina object
}
\author{Jan Oosting}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
