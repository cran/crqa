\name{checkts}
\alias{checkts}
\title{Check and trim time-series}
\description{
  Checks the length of two input time-series, and evaluates
  whether they have to be discarded given a specified
  threshold. If the difference between series is smaller
  than the threshold, the longest series between
  the two is either padded, or trimmed, to be the same length.}

\usage{
checkts(ts1, ts2, datatype, thrshd, pad)
}

\arguments{
  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{datatype}{Datatype of the time-series, either
    "numerical" or "categorical"}
  \item{thrshd}{Maximal length difference between the
    two-series used to accept them as valid.
    Series with a difference in length bigger than
    threshold are rejected as invalid.}
  \item{pad}{A logical flag indicating whether,
    in case the two-series differ by length,
    have to be padded (TRUE), with the shorter being
    extended. Or not, in which case, the longest
    series is trimmed to be of the same length
    as the shortest.}

}
\details{
  This function strictly applies when the two time-series
  series are expected to have equal length. An example of such
  case is co-registred eye-movement responses, where a
  difference in length highlights a possible error in the
  eye-tracking co-registration. The threshold is used
  to discriminate small versus large differences.
  The value to assign to the threshold should be estimated
  by looking at the difference distribution observed
  in the dataset.
  Note, when the two series are of continous numerical
  values, the padding value is obtained as the mean
  between the two series. 
}

\value{If the difference between the two-series is smaller
  than the threshold, it returns a list with two arguments:
  the first [[1]] is a matrix with two columns, i.e., the
  two time-series padded or trimmed, and number of rows equal to the
  shortest sequence. The second [[2]] is a boolean flag
  with value TRUE, indicating that the difference was
  shorter than the threshold. If instead the difference
  was bigger than threshold, it returns a list with
  two arguments: the first [[1]] is the difference observed
  between the two series, the second [[2]] is a boolean
  flag with value FALSE, indicating that the difference
  was bigger than threshold.
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}


\seealso{\code{\link{runcrqa}}}

\examples{
## generate two time-series of different length

ts1 = seq(1,30,1); ts2 = seq(1,25,1)
datatype = "continuous"
threshold = 6 ## threshold is larger than difference
pad = FALSE

res = checkts(ts1, ts2, datatype, threshold, pad)
print(res)

threshold = 4 ## threshold is smaller than difference

res = checkts(ts1, ts2, datatype, threshold)
print(res)

}

\keyword{misc}
