\name{optimizeParam}
\alias{optimizeParam}

\title{Optimal parameters value for CRQA on continuous
  time-series data}

\description{

  Iterative procedure exploring a combination of parameter
  values to obtain maximal recurrence between two
  time-series. It finds the values for the three
  parameters of radius, delay and
  embedding dimensions that optimize recurrence.
  
}
\usage{
optimizeParam(ts1, ts2, par)
}

\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{par}{A list of parameters for the optimization:
    lgM = a constant indicating maximum lag to inspect
    when calculating average mutual information between
    the two series.
    steps = a sequence of points (e.g., seq(1, 10, 1))
    used to look ahead local minima.
    cut.del = a sequence of points referring to
    the delays evaluated when mutual information
    between the two serie is estimated.} 
}

\details{

  The optimization follows a three steps process:
  
  1)  Identify a delay that accommodates both time-series
  by finding the local minimum where mutual information
  between them drops, and starts to level off.
  When one ts has a considerably longer delay indicated
  than the another, the function selects the longer delay
  of the two to ensure that new information is gained for
  both. When the delays are close to each other, the
  function computes the mean of the two delays.
  
  2)  Determine embedding dimensions by using false nearest
  neighbors and checking when it bottoms out (i.e., there
  is no gain in adding more dimensions). If the embedding
  dimension for the two ts are different the algorithm
  selects the higher embedding dimension of the two to
  make sure that both time series are sufficiently unfolded.
  
  3)  Determine radius to use for recurrence.
  Find radius that yields (1-5, percentage recurrence),
  and is in a linear scaling region of the log/log plot
  (TODO)
}

\note{As \code{optimizeParam} uses \code{crqa} to estimate
  the parameters: the additional arguments \code{normalize,
    rescale, mindiagline, minvertline, whiteline,
    recpt} should be supplied in the par list.}

\value{

  radius = r, emddim = embdim, delay = del
  
  It returns a list with the following arguments:
  \item{radius}{The optimal radius value found}
  \item{emddim}{Number of embedding dimensions}
  \item{delay}{The tau parameter for time-series delay.}
  
}
  
\references{

  Marwan, N., Carmen Romano, M., Thiel, M., and Kurths, J.
  (2007). Recurrence plots for the analysis of
  complex systems. Physics Reports, 438(5), 237-329.

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com}

\seealso{ \code{\link{crqa}}, \code{\link{wincrqa}}}

\examples{

data(crqa)

par = list(lgM =  10, steps = seq(1, 6, 1),
cut.del = seq(1, 40,1), normalize = 0,
rescale = 1, mindiagline = 2, minvertline = 2,
whiteline = FALSE, recpt = FALSE)

ans = optimizeParam(contts1, contts2, par)

print(ans)


}

\keyword{ts}
