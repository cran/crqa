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
windowdrp(x, y, step, windowsize, lagwidth, datatype, radius)

}

\arguments{

  \item{x}{First time-series}
  \item{y}{Second time-series}
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

par = list(type = 1, step =  20, windowsize =  50, lagwidth = 40,
  method = "window", datatype = "categorical", thrshd = 8)

ans = runcrqa(RDts1, RDts2, par)

print(ans)

profile = ans$profile; maxrec = ans$maxrec; maxlag = ans$maxlag

}