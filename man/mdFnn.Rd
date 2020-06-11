\name{mdFnn}
\alias{mdFnn}

\title{Find optimal embedding dimension of a multi-dimensional dataset.}

\description{Computes the percentage of false nearest
neighbors for multidimensional time series as a function 
of embedding dimension.}


\usage{ mdFnn(data, tau, maxEmb, numSamples, Rtol, Atol) }
\arguments{

 \item{data}{The matrix of data to estimate FNN.}
  \item{tau}{Time delay for embedding.}
  \item{maxEmb}{Maximum number of embedding dimensions considered}
  \item{numSamples}{Number of randomly drawn coordinates from phase-space used to estimate FNN}
  \item{Rtol}{First distance criterion for separating false neighbors}
  \item{Atol}{Second distance criterion for separating false neighbors}

}

\value{

  It returns the percentage of false neighbors for each embedding.


}

\references{Kennel, M. B., Brown, R., & Abarbanel, H. D. (1992).
Determining embedding dimension for phase-space reconstruction using
a geometrical construction. Physical review A, 45, 3403.
Wallot, S., and Moenster, D. (2018). Calculation of average mutual
information (AMI) and false-nearest neighbors (FNN) for the
estimation of embedding parameters of multidimensional time-series in
Matlab. Front. Psychol. - Quantitative Psychology and Measurement
}
\author{Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
     Dan Moenster, Aarhus University, 
     Moreno I. Coco, University of East London}

\seealso{\code{\link{mdDelay}}, \code{\link{optimizeParam}}} 

\examples{

tau = 1; maxEmb = 10; numSamples = 500; Rtol = 10; Atol = 2

data(crqa) ## load the data

handset = handmovement[1:300, ] ## take less points

mdFnn(handset, tau, maxEmb, numSamples, Rtol, Atol)

}

\keyword{ts}
