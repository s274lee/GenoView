\name{GenoView-package}
\alias{GenoView-package}
\alias{GenoView}
\docType{package}
\title{
Condensed, overlapped plotting of genomic data tracks
}
\description{
The GenoView package provides interactive visualization utilities which combine user specified and existing genomic data. 
}
\details{
\tabular{ll}{
Package: \tab GenoView\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-05\cr
License: \tab GPL-3\cr
}

The package currently focuses on visualizing mutation data with \code{\link[GenoView]{mutExonPlot}}. It provides a GUI and a text interface (see \code{\link[GenoView]{mep-Interfaces}}). Databases and variables which provide hg19 data to \code{\link[GenoView]{mutExonPlot}} are already built into \code{\link[GenoView]{mepHuman}}. 
}
\author{
Sharon Lee, Dennis Wang

Maintainer: Sharon Lee <s274lee@uwaterloo.ca>
}
\keyword{package}
\examples{
if (interactive()) {
  # Create a missense mutations dataset in TP53
  set.seed(1)
  locations = sample(7571720:7590863, size = 25)
  mut.df <- data.frame(chrom = "chr17", 
                      start = locations, 
                      end = locations, 
                      strand = "-", 
                      fill = 1:25)
  
  library(gWidgetsRGtk2)
  library(RGtk2)
  
  # hg19 GUI example
  mepHuman(dataFrame = mut.df, gui = TRUE)
}
}
