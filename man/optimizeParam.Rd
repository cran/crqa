\name{optimizeParam}
\alias{optimizeParam}

\title{Optimal parameters value for CRQA on continuous
  time-series data}

\description{

  Iterative procedure exploring a combination of parameter
  values to obtain maximal recurrence between two
  time-series. It finds the values for the three
  parameters of radius, delay and embedding 
  dimensions that optimize recurrence.
  
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
    between the two serie is estimated.
    radiusspan = a constant setting the granularity of 
    radius unit to explore, relative to the standard deviation 
    of the distance between the two series. (Larger value = smaller units) 
    radiussample = the number of equally spaced units 
    of the radius to explore. (Larger value = more units = 
    computationally more expensive) 

	} 
}

\details{

  The optimization follows a three steps process:
  
  1) Identify a delay that accommodates both time-series
  by finding the local minimum where mutual information
  between them drops, and starts to level off.
  When one ts has a considerably longer delay indicated
  than the another, the function selects the longer delay
  of the two to ensure that new information is gained for
  both. When the delays are close to each other, the
  function computes the mean of the two delays.
  
  2) Determine embedding dimensions by using false nearest
  neighbors and checking when it bottoms out (i.e., there
  is no gain in adding more dimensions). If the embedding
  dimension for the two ts are different the algorithm
  selects the higher embedding dimension of the two to
  make sure that both time series are sufficiently unfolded.
  
  3) Determine radius yielding a recurrence rate between 2-5%. 
  To do so, we first determine a starting radius that yields 
  approximately 25% recurrence rate.
  We generate a sampled sequence of equally spaced possible 
  radi from such radius till 0, using as unit for the sequence 
  step, the standard deviation of the distance matrix divided 
  by a scaling parameter (radiusspan). The larger this parameter,
  the finer the unit. The size of the sample is specified by 
  the user (radiussample). 
   
}

\note{As \code{optimizeParam} uses \code{crqa} to estimate
  the parameters: the additional arguments \code{normalize,
    rescale, mindiagline, minvertline, whiteline,
    recpt} should be supplied in the par list.
    Set up relatively large radiusspan (e.g. 100), 
    and relatively small radiussample (e.g., 10), 
    for a decent coverage of radius values.}

\value{
  
  It returns a list with the following arguments:
  \item{radius}{The optimal radius value found}
  \item{emddim}{Number of embedding dimensions}
  \item{delay}{The lag parameter.}
  
}
  
\references{

  Marwan, N., Carmen Romano, M., Thiel, M., and Kurths, J.
  (2007). Recurrence plots for the analysis of
  complex systems. Physics Reports, 438(5), 237-329.

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\seealso{ \code{\link{crqa}}, \code{\link{wincrqa}}}

\examples{

## initialize the parameters

par = list(lgM =  20, steps = seq(1, 6, 1),
cut.del = seq(1, 40,1), 
radiusspan = 100, radiussample = 40, 
normalize = 0, rescale = 1, mindiagline = 2, 
minvertline = 2, tw = 0, whiteline = FALSE, 
recpt = FALSE)

## generate two random uniform series
ts1 = runif(100)
ts2 = runif(100)

ans = optimizeParam(ts1,ts2,par)
print(ans)

## utilize leftmov, rightmov for an application to real-data
## reduce radiussample to increase computational speed

## data(crqa)
## ans = optimizeParam(leftmov, rightmov, par)
## print(ans)

}

\keyword{ts}
