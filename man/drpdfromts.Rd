\name{drpdfromts}
\alias{drpdfromts}

\title{Cross-Recurrence diagonal profile of two time-series}

\description{
  Quick method to explore the cross-recurrence diagonal
  profile of two-time series. It returns the recurrence
  observed for different delays, the maximal recurrence
  observed, and the delay at which it occurred.
}

\usage{ drpdfromts(ts1, ts2, ws, datatype, radius, 
delay, embed, rescale, normalize, mindiagline, minvertline, tw) }

\arguments{
  
  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{ws}{A constant indicating the range of delays
    (positive and negative) to explore}
  \item{datatype}{A string indicating whether the
    time-series consist of "categorical", or "continuous"
    datatype}
  \item{radius}{A threshold, cut-off, constant used
    to decide whether two points are recurrent or not.}
   \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
    if rescale = 0 (do nothing);
    if rescale = 1 (mean distance of entire matrix);
    if rescale = 2 (maximum distance of entire matrix).}
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

data(crqa)

res = drpdfromts(RDts1, RDts2, ws = 100,
datatype = "categorical", radius = 0.001, 
delay = 1, embed = 1, rescale = 0, normalize = 0, 
mindiagline = 2, minvertline = 2, tw = 0)

profile = res$profile

plot(seq(1,length(profile),1), profile, type = "l", lwd = 5,
     xaxt = "n", xlab = "Lag", ylab = "Recurrence")

}
