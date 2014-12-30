\name{mep-Interfaces}
\alias{mep-Interfaces}
\alias{mepGUI}
\alias{mepTxtInt}
\title{
Interfaces for MutExonPlot
}
\description{
Package interfaces which create, display, and save Mutations-Over-Exons plots in succession. 
}
\usage{
mepGUI(m.data, tx.db, gene.loc, seq.lens, pfam.gr = NULL, pfam.desc = NULL,
pfam.ids = NULL, ...)
mepTxtInt(m.data, tx.db, gene.loc, seq.lens, pfam.gr = NULL, pfam.desc = NULL, 
pfam.ids = NULL, ...)
}
\arguments{
  \item{m.data}{
Mutations data.frame or GRanges object which contains at minimum Chromosome, Start Position, End Position, and Strand data. It may also contain a column which uniquely identifies each mutation for legend plotting. 
}
  \item{tx.db}{
TranscriptDb object using the same reference genome as the contents of \code{m.data}
}
  \item{gene.loc}{
GRanges object containing gene data. It has HUGO gene symbols as names and has seqnames, ranges, and strand data for each entry. 
}
  \item{seq.lens}{
Seqlengths for \code{m.data}, using the same reference genome as \code{tx.db}
}
  \item{pfam.gr}{
GRanges object of protein domain coordinates and Associated identifications
}
  \item{pfam.desc}{
list of PFAM definitions with Accession numbers
}
  \item{pfam.ids}{
character vector of Associated identifications with Accession numbers as names
}
  \item{\dots}{
Input objects from other databases, datasets, vectors describing plotting options. An example of possible input objects, using the human genome, are generated in \code{\link{mepHuman}}. 
}
}
\details{
The GUI window consists of two parts, a input selection menu and a notebook widget which displays \code{m.data} and the finished plots. Once selection is complete, all required inputs are processed and passed onto \code{\link{mutExonPlot}}. The GUI displays and saves many plots at once. 

The text interface reads text entry from the user to complete input selection for plotting. The interface displays and saves one plot at a time. 
}
\author{
Sharon Lee
}
\note{
The common objects that are passed from \code{\link[GenoView]{mepHuman}} are: 
\code{tx.db}, \code{genesymbol}, \code{seqs.hg19}, \code{pfam.gr}, \code{pfam.desc}, \code{pfam.ids}, \code{int.opt}, \code{disp.opt}. 

PFAM Protein domain data might not be available for all species, and was originally obtains from UCSC. 
}

\seealso{
\code{\link[GenoView]{mepHuman}}, \code{\link[GenoView]{mutExonPlot}}
}

\examples{
\dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Create a missense mutations dataset in TP53
set.seed(1)
locations = sample(7575000:7580000, size = 25)
mut.df <- data.frame(chrom = "chr17", 
                    start = locations, 
                    end = locations, 
                    strand = "-", 
                    fill = letters[1:25])

# Generating sample data using TP53
tx.db <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqs.hg19 <- seqlengths(tx.db)
data(genesymbol, package = "biovizBase")

# Sample GRanges data with domains for P53
sample <- GRanges(seqnames = "chr17", 
                    IRanges(start = c(7576884, 7577499, 7579707), 
                            end = c(7577102, 7579403, 7579899)), 
                    strand = "-",
                    domain = c("P53_tetramer", "P53", "P53_TAD"))

# Create PFAM objects
objs <- makePFAMObjs()
desc <- objs$desc
ids <- objs$ids

int.opt <- c("Whole", "Custom")
disp.opt <- c("All", "Reduce")
dom.opt <- c("All", "Longest")

# GUI
mepGUI(m.data = mut.df, tx.db = tx.db, gene.loc = genesymbol, 
        seq.lens = seqs.hg19, pfam.desc = desc, pfam.ids = ids, 
        pfam.gr = sample, int.opt = int.opt, disp.opt = disp.opt, 
        dom.opt = dom.opt)

# Text Interface
mepTxtInt(m.data = mut.df, tx.db = tx.db, gene.loc = genesymbol, 
        seq.lens = seqs.hg19, pfam.desc = desc, pfam.ids = ids, 
        pfam.gr = sample, int.opt = int.opt, disp.opt = disp.opt, 
        dom.opt = dom.opt)
}
}