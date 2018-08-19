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
  1) it does not handle the case with more than one input,
  and  2) (m > n) matrices give the B matrix columns in a
  different order, but the d vector of indices will also
  be changed accordingly, so the set of columns is OK,
  just ordered differently
  
}

\note{For computational efficiency spdiags is actually
  computed using a Fortan implementation (\code{jspd.f}) }


\value{
  
  \item{B}{A min(m,n)-by-p matrix, usually (but not necessarily)
    full, whose columns are the diagonals of A.}
  \item{d}{A vector of length p whose integer components
    specify the diagonals in A.}
   
}
 
 
\author{John C. Nash (nashjc@uottawa.ca)}


\examples{

dta <- c(0, 5, 0, 10, 0, 0, 0, 0, 6, 0, 11, 0, 3, 0, 0,
7, 0, 12, 1, 4, 0, 0, 8, 0, 0, 2, 5, 0, 0, 9)

A1 <- matrix(dta, nrow=5, ncol=6, byrow=TRUE)

print(A1)
res1 <- spdiags(A1)
print(res1) 

}

\keyword{array}
