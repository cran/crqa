\name{theiler}
\alias{theiler}

\title{Theiler window}

\description{

  Remove recurrent points on the main diagonal, and 
  parallel diagonals, as specified by the window 
  parameter. 

}

\usage{ theiler(S, tw) }

\arguments{

  \item{S}{A binary matrix representing a recurrent plot}
  \item{tw}{The size of the theiler window.}

}


\value{

  A sparse matrix with recurrent points removed as 
  specified by the window parameter.

}

\author{Moreno I. Coco
  (moreno.cocoi@gmail.com)}

\details{The default value of the window parameter is 0.
A value of 1 would remove the points along the main diagonal.
A value of 2 would remove the points along the main diagonal,
as well as the point of the bands +|- 1 around it.}

\examples{

## build a random x matrix
r = 100; c = 100; tw = 1
S = round(matrix(runif(r*c), r, c))

ans = theiler(S, tw)

}
