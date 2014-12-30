\name{mepHuman}
\alias{mepHuman}
\title{
Main setup for Mutations Exon Plot
}
\description{
Creates all main variables using the human genome from databases and objects, passes variables and input data to interface function
}
\usage{
mepHuman(m.data, gui = FALSE)
}
\arguments{
  \item{m.data}{
Mutations data.frame or GRanges object which contains at minimum Chromosome, Start Position, End Position, and Strand data. It may also contain a column which uniquely identifies each mutation for legend plotting. 
}
  \item{gui}{
Determines whether the GUI is used. 
}
}
\details{
The mepHuman function creates the following objects which are passed to \code{mepGui} or \code{mepTxtInt}:

\code{tx.db} UCSC hg19 known Gene database from \code{\link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}}

\code{genesymbol} GRanges object describing human genes from \code{\link[biovizBase]{genesymbol}}

\code{seqs.hg19} seqlengths of hg19 object

\code{pfam.gr}, \code{pfam.desc}, \code{pfam.ids} objects containing PFAM domain information

\code{int.opt}, \code{disp.opt} character vectors which describe the plotting input choices
}
\value{
A mutExonPlot using the hg19 platform
}
\author{
Sharon Lee
}
\seealso{
\code{\link{mutExonPlot}}, \code{\link{mep-Interfaces}}
}
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
  
  # GUI example
  mepHuman(m.data = mut.df, gui = TRUE)
}

if (interactive()) {
  # Text interface example
  mepHuman(m.data = mut.df)
}
}
