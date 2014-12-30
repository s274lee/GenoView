\name{PFAMIDE}
\alias{PFAMIDE}
\title{
Translates domain transcripts' identifications to descriptions
}
\description{
Maps the domain values of a GRanges object from Associated identifications to PFAM definitions
}
\usage{
PFAMIDE(transcripts, desc, ids)
}
\arguments{
  \item{transcripts}{
GRanges object with \code{domain} column in its values containing PFAM Associated identifications as characters
}
  \item{desc}{
list of PFAM definitions with Accession numbers
}
  \item{ids}{
character vector of Associated identifications with Accession numbers as names
}
}
\value{
GRanges object with mapped \code{domain} column in its values
}
\author{
Sharon Lee
}
\note{
There are some Accession numbers from \code{ids} that do not map to those in \code{desc}. The original Associated identification is retained, rather than producing \code{NA}. 
}

\seealso{
\code{\link[GenoView]{makePFAMObjs}}
}
\examples{
\dontrun{
# Sample data with domains for P53
sample <- GRanges(seqnames = "chr17", 
                    IRanges(start = c(7576884, 7577499, 7579707), 
                            end = c(7577102, 7579403, 7579899)), 
                    strand = "-",
                    domain = c("P53_tetramer", "P53", "P53_TAD"))

# Create the mapping objects
pfam.objs <- makePFAMObjs()

PFAMIDE(transcripts = sample, desc = pfam.objs$desc, ids = pfam.objs$ids)
}
}
