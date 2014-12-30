\name{colSel}
\alias{colSel}
\title{
Column name selection
}
\description{
Selects a column name by user input.
}
\usage{
colSel(col.crit, dataFrame)
}
\arguments{
  \item{col.crit}{
Description of type of data in the desired column.
}
  \item{dataFrame}{
Matrix or data.frame with column names.
}
}
\value{
Character which is a column name in \code{dataFrame}
}

\author{
Sharon Lee
}

\seealso{
\code{\link[GenoView]{printOption}}
}
\examples{
if(interactive()) {
  colSel("Gear", dataFrame = mtcars)
}
}
