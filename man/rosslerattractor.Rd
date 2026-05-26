\name{rosslerattractor}
\alias{rosslerattractor}
\title{Simulate a Roessler attractor}
\description{
Companion to \code{\link{lorenzattractor}}. Generates a numerical
trajectory of the Roessler system via plain Euler integration; mainly
used to drive the package's benchmark / validation scripts and to
provide a second canonical chaotic test system.
}
\usage{
rosslerattractor(numsteps, dt = 0.05, a = 0.25, b = 0.25, c = 4)
}
\arguments{
  \item{numsteps}{Number of integration steps to produce.}
  \item{dt}{Integration time-step.}
  \item{a, b, c}{Roessler system parameters.
    The defaults reproduce the parameter set used in
    Marwan & Kraemer (2023, EPJST 232:5-27, Fig. 3A and
    Appendix A.3) for benchmarking RQA packages.}
}
\details{
Integrates the system
\deqn{\dot x = -y - z, \quad \dot y = x + a y, \quad \dot z = b + z (x - c)}
with a small random initial condition (\code{rnorm(3)}). Discard the
first 1000 points as transient before using the trajectory for RQA.
}
\value{
A numeric matrix with \code{numsteps} rows and three columns
(\code{x}, \code{y}, \code{z}).
}

\references{
Roessler, O. E. (1976). An equation for continuous chaos. Phys. Lett. A 57(5):397-398.

Marwan, N., Kraemer, K. H. (2023). Trends in recurrence analysis of dynamical systems. Eur. Phys. J. Spec. Top. 232:5-27.
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\seealso{\code{\link{lorenzattractor}}}

\examples{
set.seed(42)
ross <- rosslerattractor(2000, dt = 0.05)
ts   <- ross[1001:2000, 1]   # drop transient, keep x-component
plot(ts, type = "l", main = "Roessler x-component")
}

\keyword{datagen}
