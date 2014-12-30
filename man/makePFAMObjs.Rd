\name{makePFAMObjs}
\alias{makePFAMObjs}
\title{
Make PFAM mapping objects
}
\description{
Produces custom objects from PFAM.db, targets PFAM mapping in \code{\link[GenoView]{mutExonPlot}}
}
\usage{
makePFAMObjs()
}
\value{
  \item{desc}{
list of PFAM definitions with Accession numbers
}
  \item{ids}{
character vector of Associated identifierss with Accession numbers as names
}
}
\author{
Sharon Lee
}
\examples{
makePFAMObjs()
}