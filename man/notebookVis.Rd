\name{notebookVis}
\alias{notebookVis}
\title{
Visualize Grob objects in gnotebook widget
}
\description{
Create a new page in the gnotebook widget to display the plot with a unique label
}
\usage{
notebookVis(nb, plt, type, suff)
}
\arguments{
  \item{nb}{
Visible gnotebook widget
}
  \item{plt}{
Grob object to plot
}
  \item{type}{
Description of \code{plt} in tab label
}
  \item{suff}{
Incremental numeric to follow \code{type} in tab label
}
}
\value{
\item{plt}{
The same input \code{plt} which has been visualized
}
\item{plt.lab}{
The tab label for the notebook page displaying \code{plt}
}
}
\author{
Sharon Lee
}
\examples{
if (interactive()) {
library(grid)
library(gWidgets)

# Create a grob object
l <- linesGrob()

# Create the notebook widget
win <- gwindow()
nb <- gnotebook(cont = win)
notebookVis(nb = nb, plt = l, type = "Example", suff = 1)
}
}