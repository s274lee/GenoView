\name{printOption}
\alias{printOption}
\title{
Print numbered options
}
\description{
Prints a numbered list as a guide for selecting an option or an index. 
}
\usage{
printOption(opt.vec, inst = TRUE, desc = "an option", 
match.opt = FALSE, match.ind = FALSE)
}
\arguments{
  \item{opt.vec}{
A vector of printable objects, preferably character, which fits the descriptor \code{desc}.
}
  \item{inst}{
Controls instruction printing using \code{desc}.
}
  \item{desc}{
A char which indicates the type of object in \code{opt.vec}.
}
  \item{match.opt}{
Chooses the selected element from \code{opt.vec}.
}
  \item{match.ind}{
Chooses the index of the selected element from \code{opt.vec}.
}
}
\value{
The selected object from \code{opt.vec} if \code{match.opt = TRUE}.
The integer index if \code{match.ind = TRUE}.
}
\author{
Sharon Lee
}
\examples{
colours <- c("Red", "Blue", "Green", "Yellow")
printOption(colours, desc = "a colour") ## no return value

if (interactive()){
printOption(colours, desc = "a colour", match.ind = TRUE) # returns integer

printOption(colours, desc = "a colour", match.opt = TRUE) # returns character
}
}