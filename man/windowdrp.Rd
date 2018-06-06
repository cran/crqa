\name{windowdrp}
\alias{windowdrp}

\title{Window Cross-Recurrence Profile}

\description{
  
  Cross-recurrence is calculated in overlapping
  windows of a specified size for a number of
  delays smaller than the size of the window.
  In every window, the recurrence value for the
  different delays is calculated.
  A mean is then taken across the delays to obtain
  a recurrence value in that particular window.

}

\usage{
windowdrp(ts1, ts2, step, windowsize, lagwidth, datatype, 
radius, delay, embed, rescale, normalize, mindiagline, minvertline,
tw)

}

\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{step}{Interval by which the window is
    moved.}
  \item{windowsize}{The size of the window}  
  \item{lagwidth}{The number of delays to be
    considered}
  \item{datatype}{A string indicating whether
    the time-series is 'numeric' or 'categorical'}
  \item{radius}{For numeric time-series, the
    cutoff distance to accept or reject
    two-points as recurrent}
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
  
  It returns a list of arguments where:
  \item{profile}{Time-course cross-recurrence profile}
  \item{maxrec}{Maximal recurrence observed along the time-series}
  \item{maxlag}{The point in time where maximal
    recurrence is observed}

}

\references{

  Boker, S. M., Rotondo, J. L., Xu, M.,
  & King, K. (2002).
  Windowed cross-correlation and peak picking
  for the analysis of variability in the
  association between behavioral time series.
  Psychological Methods, 7(3), 338.

}
\author{Moreno I. Coco (moreno.cocoi@gmail.com)
  and Rick Dale (rdale@ucmerced.edu)}


\seealso{\code{\link{drpdfromts}}}

\examples{

data(crqa)

ts1 = as.numeric(as.matrix(RDts1))
ts2 = as.numeric(as.matrix(RDts2))

par = list(type = 1, step =  20, windowsize =  50, lagwidth = 40,
           method = "window", datatype = "categorical", thrshd = 8,
           radius = 0.001, delay = 1, embed = 1, rescale = 0, normalize = 0,
           mindiagline = 2, minvertline = 2)

ans = runcrqa(ts1, ts2, par)

print(ans)

profile = ans$profile; maxrec = ans$maxrec; maxlag = ans$maxlag

}