\name{sample.mut.df}
\alias{sample.mut.df}
\docType{data}
\title{
Sample hg19 TP53 missense mutations dataset
}
\description{
Data.frame containing hg19 human genome TP53 missense mutations, to be used in examples and vignette
}
\usage{data(sample.mut.df)}
\format{
  A data frame with 10 observations on the following 7 variables.
  \describe{
    \item{\code{ids}}{a factor with sample ID}
    \item{\code{seqnames}}{a factor with chromosome locations}
	\item{\code{start}}{a numeric, integer vector}
    \item{\code{end}}{a numeric, integer vector}
    \item{\code{width}}{a numeric, integer vector}
    \item{\code{strand}}{a factor with levels \code{-} \code{+}}
    \item{\code{midpoint}}{a numeric vector}
  }
}

\examples{
\dontrun{
data(sample.mut.df)
}
}
\keyword{datasets}
