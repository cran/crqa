\name{optimizeParam}
\alias{optimizeParam}

\title{Estimate optimal delay, embedding dimension and radius
for continuous time-series data}

\description{

  Iterative procedure to examine the values
  of delay, embedding dimension and radius 
  to compute recurrence plots of one, two, 
  or more time-series.
  
}
\usage{
optimizeParam(ts1, ts2, par, min.rec, max.rec)
}

\arguments{
  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{par}{A list of parameters needed for the optimization, 
  refer to the Details section.}
  \item{min.rec}{The minimum value of recurrence accepted.
    Default = 2}
  \item{max.rec}{The maximum value of recurrence accepted.
    Default = 5}
}

\details{

  The optimization can be applied both to uni-dimensional 
  time-series (method = crqa), or multi-dimensional (method = mdcrqa)
  
  The procedure is identical in both cases:
  
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
  the finer the unit.
  
  For uni-dimensional time-series, the user has to decide how to 
  choose the value of average mutual information (i.e., \code{typeami = mindip}, 
  the lag at which minimal information is observed, or \code{typeami = maxlag}, 
  the maximum lag at which minimal information is observed) and the relative
  percentage of information gained in FNN, relative to the first embedding 
  dimension, when higher embeddings are considered (i.e., \code{fnnpercent}). 
  Then, as \code{crqa} is integrated in the \code{optimizeParam} to estimate the
  radius, most of the arguments are the same (e.g., \code{mindiagline} 
  or \code{tw}).
  
  For multidimensional series, the user needs to specify the right RQA method 
  (i.e., \code{method = "mdcrqa"}). Then, for the estimation of the delay via
  AMI: (1) \code{nbins} the number of bins to compute the two-dimensional 
  histogram of the original and delayed time series and (2) the \code{criterion} 
  to select the delay (\code{firstBelow} to use the lowest delay at which the AMI
  function drops below the value set by the \code{threshold} argument, and 
  \code{localMin} to use the position of the first local AMI minimum). 
  The estimation of the embedding dimensions instead needs the following arguments:
  (1) \code{maxEmb}, which is the maximum number of embedding dimensions considered,
  (2) \code{noSamples}, which is the number of randomly drawn coordinates 
  from phase-space used to estimate the percentage of false-nearest neighbors,
  (3) \code{Rtol}, which is the first distance criterion for separating false
  neighbors, and (4) \code{Atol}, which is the second distance criterion for
  separating false neighbors. The radius is estimated as before.
   
}

\note{As \code{optimizeParam} uses \code{crqa} to estimate
  the parameters: the additional arguments normalize,
    rescale, mindiagline, minvertline, whiteline,
    recpt should be supplied in the par list.
    Set up relatively large radiusspan (e.g. 100), 
    for a decent coverage of radius values.
  }

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

\author{Moreno I. Coco (moreno.cocoi@gmail.com),
James A. Dixon (james.dixon@uconn.edu)
Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
Dan Moenster, Aarhus University
}

\seealso{ \code{\link{crqa}}, \code{\link{wincrqa}}}

\examples{

data(crqa) ## load the data

handset = handmovement[1:300, ] ## take less points

P1 = cbind(handset$P1_TT_d, handset$P1_TT_n) 
P2 = cbind(handset$P2_TT_d, handset$P2_TT_n)

par = list(method = "mdcrqa", metric = "euclidean", maxlag =  20, 
           radiusspan = 100, radiussample = 40, normalize = 0, 
           rescale = 4, mindiagline = 10, minvertline = 10, tw = 0, 
           whiteline = FALSE, recpt = FALSE, side = "both", 
           datatype = "continuous", fnnpercent  = NA,  
           typeami = NA, nbins  = 50, criterion = "firstBelow",
           threshold = 1, maxEmb = 20, numSamples = 500, 
           Rtol = 10, Atol = 2)

results = optimizeParam(P1, P2, par, min.rec = 2, max.rec = 5)
print(unlist(results))
           

}

\keyword{ts}
