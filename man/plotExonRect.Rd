\name{plotExonRect}
\alias{plotExonRect}
\title{
Plot Rectangular Exon Components
}
\description{
Adds colored rectangular blocks, denoting exon components such as coding sequences or untranslated regions, to a ggplot. Provides filtering to highlight a subset of the input data in the plot. 
}
\usage{
plotExonRect(plt, component = c("cds", "utr"), comp.df,
col.name = "type", filt.by.comp = FALSE, size = c(0.4, 0.15), 
subset.df, comp.col = "dark grey", 
sub.col = RColorBrewer::brewer.pal(3, "Set1")[2], 
aesCst)
}
\arguments{
  \item{plt}{
A ggplot object
}
  \item{component}{
Exon component string, such as "cds", "utr"
}
  \item{comp.df}{
Exon components data.frame, contains a \code{col.name} column which contains exon component strings and columns for start, end, and stepping (starting y position in plot)
}
  \item{col.name}{
Column name used to filter \code{comp.df} before plotting
}
  \item{filt.by.comp}{
Determines if \code{comp.df} should be filtered for \code{component} contents before plotting
}
  \item{size}{
Half the height of the rectangular blocks
}
  \item{subset.df}{
Secondary data.frame, contains a \code{col.name} column which contains exon component strings and columns for start, end, and stepping (starting y position in plot). This data is plotted in a highlighting colour. 
}
  \item{comp.col}{
Colour of data in \code{code.df}
}
  \item{sub.col}{
Highlighting colour of data in \code{subset.df}
}
  \item{aesCst}{
Custom aes function which can handle constants in variables
}
}
\value{
ggplot object with additional geom_rect() layers
}
\author{
Sharon Lee
}
\examples{
library(ggplot2)

set.seed(3)
N = 30
# Create sample base exon components
starts = sample(1:10000, size = N, replace = TRUE)
sample.df <- data.frame(start = starts, 
                        end = starts + sample(100:500, size = N), 
                        stepping = 1, 
                        type = sample(c("cds", "utr", "gap"), 
                                        size = N, replace = TRUE))

# Create sample highlighted exon components
N = 5
starts = sample(1:10000, size = N, replace = TRUE)
sample.sub <- data.frame(start = starts, 
                        end = starts + sample(100:500, size = N), 
                        stepping = 1, 
                        type = sample(c("cds", "utr", "gap"), 
                                        size = N, replace = TRUE))

# Using an aes function which evaluates arguments locally without parsing
aesCst <- function(...) {
    structure(list(...),  class = "uneval")
}

step1 <-ggplot() + 
        geom_abline(intercept = 1, slope = 0, 
                    aes(xmin = -1, xmax = 1, ymin = 0.5, ymax = 1.5), 
                    color = "dark grey")
                    
step2 <- plotExonRect(plt = step1, component = "cds", comp.df = sample.df, 
                        filt.by.comp = TRUE, size = 0.4, subset.df = sample.sub, 
                        aesCst = aesCst)
plotExonRect(plt = step2, component = "utr", comp.df = sample.df, 
                        filt.by.comp = TRUE, size = 0.15, subset.df = sample.sub, 
                        aesCst = aesCst)
}
