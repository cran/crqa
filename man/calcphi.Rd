\name{calcphi}
\alias{calcphi}

\title{Extract the phi-coefficient observed between the 
two time-series on a specific state k.}

\description{Phi-coefficient is the recurrence observed between the two
  time-series on a specific state k. The phi(k) coefficient
  increases with the frequency of matching recurrence on the same
  state (k ; k) and away from this state (not k ; not k) between
  the two time-series. On the other hand, phi(k) decreases with the
  frequency of mismatching objects (k; not-k, and vice versa).}

\usage{ calcphi(t1, t2, ws, k) }
\arguments{
  \item{t1}{First time-series}
  \item{t2}{Second time-series}
  \item{ws}{Number of delays (+/-) considered}
  \item{k}{The categorical state on which phi is calculated}
}

\value{

  It returns the recurrence phi-coefficient profile for
  state k for all delays considered

}

\references{Dale, R., Warlaumont, A. S. and Richardson, D. C. (2011).
  Nominal cross recurrence as a generalized lag sequential analysis for
  behavioral streams. International Journal of Bifurcation and Chaos,
  21, 1153-1161. (special issue on recurrence)}
\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\seealso{\code{\link{CTcrqa}}, \code{\link{simts}}} 

\examples{

## simulate two dichotomous series
tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 100)
ts1 = tS[1,]; ts2 = tS[2,]

## k = 1, as series are dichotomous
## check data(crqa) for alternative data (RDts1, RDts2)

k = 1; ws = 40
res = calcphi(ts1, ts2, ws, k)

}

\keyword{ts}
