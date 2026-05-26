\name{spdiags}
\alias{spdiags}

\title{Extract diagonal matrices}

\description{
  Extracts all nonzero diagonals from the m-by-n matrix A.
  B is a min(m,n)-by-p matrix whose columns are the
  p nonzero diagonals of A.

}
\usage{ spdiags(A) }

\arguments{

  \item{A}{An m-by-n matrix with nonzero elements located on p diagonals.}
  
}

\details{

  Compared to the original Matlab implementation:
  1) it does not handle the case with more than one input.
  2) Columns in \code{B} are always ordered by ascending diagonal index
  (\code{d} is sorted increasingly) for all matrix shapes.
  3) An all-zero input returns an empty \code{B} (0 columns) consistent
  with the contract of extracting \emph{nonzero} diagonals.

}

\note{Implemented in pure R using vectorized triplet indexing.
  No compiled code is required.}


\value{
  
  \item{B}{A min(m,n)-by-p matrix, usually (but not necessarily)
    full, whose columns are the diagonals of A.}
  \item{d}{A vector of length p whose integer components
    specify the diagonals in A.}
   
}
 
 
\author{Moreno I. Coco (moreno.cocoi@gmail.com)}


\examples{

dta <- c(0, 5, 0, 10, 0, 0, 0, 0, 6, 0, 11, 0, 3, 0, 0,
7, 0, 12, 1, 4, 0, 0, 8, 0, 0, 2, 5, 0, 0, 9)

A1 <- matrix(dta, nrow=5, ncol=6, byrow=TRUE)

print(A1)
res1 <- spdiags(A1)
print(res1) 

}

\keyword{array}
