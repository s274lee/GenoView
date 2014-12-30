\name{mutExonPlot}
\alias{mutExonPlot}
\title{
Mutations-Over-Exons Plot
}
\description{
Superimpose user specified mutation data on the exons of a gene, with domain and legend annotations
}
\usage{
mutExonPlot(mut.gr, exon.int, plot.int = exon.int, 
disp.track = 2, p.height = 1/4, d.height = 0, l.height = 1/2, 
plt.title = "mutExonPlot", id.col = "Fill", tx.db = NULL, pfam.gr = NULL, 
pfam.desc = NULL, pfam.ids = NULL, gr = NULL, ...)
}
\arguments{
  \item{mut.gr}{
GRanges object containing mutation data, with fill values for legend if applicable
}
  \item{exon.int}{
GRanges object of length 1 of the genomic interval in which to count exons
}
  \item{plot.int}{
GRanges object of length 1 of the plotting interval in genomic coordinates
}
  \item{disp.track}{
Selects exon plotting track: 1 for full transcripts, 2 for reduced track
}
  \item{p.height}{
Numeric height allocation of Mutations-Over-Exons section in final grob 
}
  \item{d.height}{
Numeric height allocation of domain annotation section in final grob
}
  \item{l.height}{
Numeric height allocation of legend section in final grob
}
  \item{plt.title}{
Plot title
}
  \item{id.col}{
Column containing identifiers in the values of \code{mut.gr}. It is only relevant in plotting a legend. 
}
  \item{tx.db}{
TrancscriptDb object
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
  \item{gr}{
GRanges object containing genomic data (cds, utr, gap, intron, exon) over \code{exon.int}
}
  \item{\dots}{
Objects passed from the interface function if applicable
}
}
\value{
Some elements may not be returned, depending on the chosen plot components.
\item{mep}{
Mutations-Over-Exons final plot grob
}
\item{plot}{
Mutations-Over-Exons component grob
}
\item{domain}{
Domain component grob
}
\item{legend}{
Legend component grob
}
\item{show.domain}{
show.domain boolean
}
\item{gr}{
GRanges object containing genomic data (cds, utr, gap, intron, exon) over \code{exon.int}
}
\item{exon.int}{
GRanges object of length 1 which specifies the gene's genomic coordinates 
}

If plotting fails, \code{FALSE} is returned. 
}
\author{
Sharon Lee
}
\note{
\code{plot.int} should be smaller than or equal to \code{exon.int}. 

When the reduced track option and domain plotting are selected in \code{disp.track} and \code{d.height}, the longest protein domain in the plotting interval is added the mutations over exon ploting section, rather than its own domain section. 
}

\examples{
\dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biovizBase)

# Determine the published genomic datasets
tx.db <- TxDb.Hsapiens.UCSC.hg19.knownGene
data(genesymbol, package = "biovizBase")

# Create a mutation dataset in TP53
exon.int <- genesymbol["TP53"]
set.seed(2)
N = 10
mut.gr <- GRanges(seqnames = "chr17", 
                    IRanges(start = sample(7575000:7580000, 
                            size = N, replace = TRUE), 
                            width = 1), 
                    strand = "-")

mutExonPlot(mut.gr = mut.gr, exon.int = exon.int, l.height = 0, tx.db = tx.db)
}
}