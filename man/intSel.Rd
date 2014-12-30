\name{intSel}
\alias{intSel}
\title{
GRanges interval selection
}
\description{
Truncate the range of a GRanges object
}
\usage{
intSel(gr, new.start = NULL, new.end = NULL)
}
\arguments{
  \item{gr}{
GRanges object of length 1
}
  \item{new.start}{
New start location 
}
  \item{new.end}{
New end location
}
}
\details{
The start and end locations, input or user selected, create an interval which is smaller or equal to the original range of \code{gr}
}
\value{
A Granges object of length 1
}
\author{
Sharon Lee
}
\examples{
a <- GRanges(seqnames = "chr1", 
                ranges = IRanges(start = 1000, end = 5000), 
                strand = "-")
if (interactive()) {
  intSel(a)
}

b <- GRanges(seqnames = "chr2", 
                ranges = IRanges(start = 10000, end = 50000), 
                strand = "-")
intSel(gr = b, new.start = 12345, new.end = 49999)
}