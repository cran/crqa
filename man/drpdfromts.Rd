\name{drpdfromts}
\alias{drpdfromts}

\title{Cross-Recurrence diagonal profile of two time-series}

\description{
  Quick method to explore the cross-recurrence diagonal
  profile of two-time series. It returns the recurrence
  observed for different delays, the maximal recurrence
  observed, and the delay at which it occurred.
}

\usage{ drpdfromts(t1, t2, ws, datatype, radius) }

\arguments{
  
  \item{t1}{First time-series}
  \item{t2}{Second time-series}
  \item{ws}{A constant indicating the range of delays
    (positive and negative) to explore}
  \item{datatype}{A string indicating whether the
    time-series consist of "categorical", or "continuous"
    datatype}
  \item{radius}{A threshold, cut-off, constant used
    to decide whether two points are recurrent or not.}
  
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
datatype = "categorical", radius = 0.000001)

profile = res$profile

plot(seq(1,length(profile),1),profile,type = "l", lwd = 5,
     xaxt = "n", xlab = "Delays", ylab = "Recurrence")

}
