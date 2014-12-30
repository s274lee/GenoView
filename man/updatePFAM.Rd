\name{updatePFAM}
\alias{updatePFAM}
\title{
Process raw UCSC PFAM table
}
\description{
Truncates and maps the contents of the raw hg19 UCSC PFAM table
}
\usage{
updatePFAM(raw.file = "ucscGenePfam.txt", search.dir = getwd(), 
file.name = "", ID2DE = FALSE)
}
\arguments{
  \item{raw.file}{
Filename of raw PFAM table
}
  \item{search.dir}{
Directory in which the raw.file should be located
}
  \item{file.name}{
Filename of the processed table
}
  \item{ID2DE}{
Controls mapping from PFAM identification for Accession number to its definition
}
}
\details{
The protein domain table can be obtained by downloading from the hg19 annotation database of UCSC (see default \code{raw.file}), or extracting data using the UCSC Table Browser Tool. 
}
\value{
Produces a truncated tab-delimited text file named \code{file.name} in \code{search.dir}. 

Returns the processed file name if a file was created, FALSE otherwise. 
}
\references{
Karolchik D, et al. The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004;32:D493-D496.
}
\author{
Sharon Lee
}
\examples{
updatePFAM()
}