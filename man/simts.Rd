\name{simts}
\alias{simts}

\title{Simulate dichotomous binary time-series}

\description{

  A simple algorithm for producing a time-series that drives
  a second time-series (1 for event occurrence; 0 otherwise)
  using parameters, which change independent and conditional
  probability of an event to occur.

  
}

\usage{ simts(BL1, BL2, BLR1, BLR2, BL2C1, tsL) }

\arguments{

  \item{BL1}{Base event rate of the first time-series}
  \item{BL2}{Base event rate of the second time-series}
  \item{BLR1}{Rate of repetition in the first series}
  \item{BLR2}{Rate of repetition in the second series}
  \item{BL2C1}{Conditional probability of repetition.}
  \item{tsL}{Length of the simulated time-series}

}


\value{

  A matrix with two-rows, where the first row is
  the 'driving` time-series and the second row
  is the second time-series. The columns are
  the number of simulated points as selected
  by the argument \code{tsL}.

}

\author{Rick Dale and Moreno I. Coco
  (moreno.cocoi@gmail.com)}


\examples{

## set up parameters

BL1 = .08; BL2 = .05; BLR1 = .5; BLR2 = .5;
BL2C1 = .33; tsL = 100

ts = simts(BL1, BL2, BLR1, BLR2, BL2C1, tsL)

}

\keyword{ts}
