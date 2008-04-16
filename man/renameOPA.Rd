\name{renameOPA}
\alias{renameOPA}
\title{Change the linkage panel in a dataset}
\description{
  Change the linkage panel in a dataset
}
\usage{
renameOPA(snpdata, newOPA)
}
\arguments{
  \item{snpdata}{SnpSetIllumina object}
  \item{newOPA}{character, new linkage panel}
}
\details{
  In order to combine different versions of the linkage panels, this function makes it 
  possible to map the equivalent SNPs in both datasets.
}
\value{
  SnpSetIllumina object
}
\author{Jan Oosting}
\keyword{manip}
