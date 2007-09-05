\name{alterCN}
\alias{alterCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ alterCN }
\description{
  Changes one of the levels of a cn.sum data structure
}
\usage{
alterCN(cn.sum, opa, value, updown)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cn.sum}{cn.sum structure to change}
  \item{opa}{opa panel within the structure}
  \item{value}{the predicted value to change}
  \item{updown}{the value has a higher (TRUE) or lower (FALSE) cn value}
}
\details{
  The state in the cn.sum structure that has a predicted value of \code{value} 
  will have it's associated associated inferred copy number increased ( \code{updown}
  is \code{TRUE}) or decreased (\code{updown} is \code{FALSE}).
  The function makes sure that the copynumber values within a OPA panel have the 
  same order as the predicted values. 
  
}
\value{
  a new cn.sum data structure
}
\author{Jan Oosting}
\note{ 
}
\seealso{ \link{interactiveCNselect}, \link{createCNSummary}, \link{setRealCN} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
