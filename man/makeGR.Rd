\name{makeGR}
\alias{makeGR}
\title{
Make GRanges object from data.frame
}
\description{
Extracts GRanges data from relevant columns from a data.frame, performs seqlengths adjustments and complete case filtering
}
\usage{
makeGR(dataFrame, chr.col, s.p.col, e.p.col, str.col, 
id.col = NULL, plot.int, des.seqs, show.legend = FALSE, l.height = 0)
}
\arguments{
  \item{dataFrame}{
Data.frame which contains at least Chromosome, Start Position, End Position, and Strand data. It may also contain a column which uniquely identifies each entry for legend plotting. 
}
  \item{chr.col}{
Chromosome data column name
}
  \item{s.p.col}{
Start position data column name
}
  \item{e.p.col}{
End position data column name
}
  \item{str.col}{
Strand data column name
}
  \item{id.col}{
Fill data column name
}
  \item{plot.int}{
GRanges object of length 1 of the plotting interval in genomic coordinates
}
  \item{des.seqs}{
Desired seqlengths appropriate for \code{dataFrame}
}
  \item{show.legend}{
Determines whether legend information is kept based on boolean
}
  \item{l.height}{
Determines whether legend information is kept based on height allocation
}
}
\value{
GRanges object containing filtered genomic data
}
\author{
Sharon Lee
}
\examples{
library(biovizBase)

# Create a missense mutations dataset in TP53
set.seed(1)
N = 25
locations = sample(7571720:7590863, size = N)
mut.df <- data.frame(chrom = "chr17", 
                    start = locations, 
                    end = locations, 
                    str = "-", 
                    fill = 1:N)

# Create a truncated plotting interval
plot.int = GRanges(seqnames = "chr17", 
                    IRanges(start = 7576000, end = 7579000), 
                    strand = "-")

# Obtain hg19 seqlengths data
data(hg19Ideogram, package = "biovizBase")
seqs.hg19 <- seqlengths(hg19Ideogram)

makeGR(dataFrame = mut.df, chr.col = "chrom", s.p.col = "start", 
        e.p.col = "end", str.col = "str", id.col = "fill", 
        plot.int = plot.int, show.legend = TRUE, 
        des.seqs = seqs.hg19)
}