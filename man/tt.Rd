\name{tt}
\alias{tt}

\title{Trapping-Time}

\description{

  Extract vertical lines from a recurrence plot
  on which it calculates laminarity (the percentage of
  recurrence points which form vertical lines),
  and trapping-time (the mean length of vertical lines).
  
}

\usage{ tt(x, minvertline, whiteline) }

\arguments{

  \item{x}{A binary matrix representing a recurrent plot}
  \item{minvertline}{A minimum vertical length of
    recurrent points.}
  \item{whiteline}{A logical flag to calculate (TRUE)
    or not (FALSE) empty vertical lines}

}


\value{

  A list with four variables: TT (The average length of vertical
  line structures), lam (proportion of recurrent points forming
  vertical line structures), tw (vertically consecutive white
  points|lines) and tb (vertically consecutive black points|lines).
  If whiteline = FALSE, tw = NA.

}

\author{Moreno I. Coco
  (moreno.cocoi@gmail.com)}

\details{This function is based on original MATLAB code
  written by Norbert Marwan, available in \code{crptoolbox}}

\examples{

## build a random x matrix
r = 100; c = 100
x = round(matrix(runif(r*c), r, c))
whiteline = TRUE
minvertline = 2

ans = tt(x, minvertline, whiteline)

}
