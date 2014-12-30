\name{pfam.df}
\alias{pfam.df}
\docType{data}
\title{
Human protein domain dataset
}
\description{
Processed dataset containing PFAM protein domain locations and identifications
}
\usage{data(pfam.df)}
\format{
  A data frame with 343252 observations on the following 5 variables.
  \describe{
    \item{\code{PFAMID}}{a factor with PFAM Associated identifications as levels}
	\item{\code{chr}}{a factor with chromosomes as levels}
	\item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{strand}}{a factor with levels \code{-} \code{+}}
  }
}
\source{
UCSC Table Browser hg19 known to Pfam table (default file name "ucscGenePfam.txt"") retrieved 2014/12/29
}
\references{
Karolchik D, et al. The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004;32:D493-D496.
}
\examples{
data(pfam.df)
}
\keyword{datasets}
