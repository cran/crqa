\name{drpfromts}
\alias{drpfromts}

\title{Diagonal recurrence profile}

\description{
  Method to explore the diagonal profile of 
  the recurrence plot (Auto, Cross, or Multi-dimensional).
  It returns the recurrence for different delays, 
  the maximal recurrence observed and the delay at which it occurred.
}

\usage{ drpfromts(ts1, ts2, windowsize, radius,
delay, embed, rescale, normalize, mindiagline, minvertline, tw,
whiteline, recpt, side, method, metric, datatype) }

\arguments{
  
  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{windowsize}{A constant indicating the range of delays
    (positive and negative) to explore}
  \item{radius}{A threshold, cut-off, constant used
    to decide whether two points are recurrent or not.}
   \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
  if rescale = 0 (do nothing);
    if rescale = 1 (mean distance of entire matrix);
    if rescale = 2 (maximum distance of entire matrix).
    if rescale = 3 (minimum distance of entire matrix).
    if rescale = 4 (euclidean distance of entire matrix).}
  \item{normalize}{Normalize the time-series;
    if normalize = 0 (do nothing);
    if normalize = 1 (Unit interval);
    if normalize =  2 (z-score).}
  \item{mindiagline}{A minimum diagonal length of
    recurrent points. Usually set to 2, as it takes
    a minimum of two points to define any line.}
  \item{minvertline}{A minimum vertical length of
    recurrent points.}
  \item{tw}{The Theiler window parameter}
  \item{whiteline}{A logical flag to calculate (TRUE)
    or not (FALSE) empty vertical lines.}
  \item{recpt}{A logical flag indicating whether
    measures of cross-recurrence are calculated directly
    from a recurrent plot (TRUE) or not (FALSE).}
  \item{side}{A string indicating whether recurrence measures
    should be calculated in the 'upper' triangle of the RP
    'lower' triangle of the matrix, or 'both'.
    LOC is automatically excluded for 'upper' and 'lower'.}
  \item{method}{A string to indicate the type of recurrence
  analysis to peform. There are three options:
  rqa (autorecurrence); crqa(cross-recurrence); 
  mdcrqa(multidimensional recurrence). Default value is crqa}
  \item{metric}{A string to indicate the type of 
  distance metric used, default is euclidean but 
  see help rdist() to list all other possible metrics.} 
  \item{datatype}{a string (continuous or categorical) 
  to indicate whether the nature of the data type}
  }

\value{
  A list with the following arguments:
  \item{profile}{A vector of recurrence (ranging from 0,1)
    with length equal to the number of delays explored}
  \item{maxrec}{Maximal recurrence observed between the
    two-series}
  \item{maxlag}{Delay at which maximal recurrence is
    observed}
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}


\seealso{\code{\link{windowdrp}}}

\examples{

# use the available data
data(crqa) 

listener = eyemovement$listener
narrator = eyemovement$narrator

res = drpfromts(narrator, listener, windowsize = 100,
                 radius = 0.001, delay = 1, embed = 1, rescale = 0,
                 normalize = 0, mindiagline = 2, minvertline = 2,
                 tw = 0, whiteline = FALSE, recpt = FALSE, 
                 side = 'both', method = 'crqa', 
                 metric = 'euclidean', datatype = 'continuous')

 profile = res$profile

 plot(seq(1,length(profile),1), profile, type = "l", lwd = 5,
     xaxt = "n", xlab = "Lag", ylab = "Recurrence")

}
