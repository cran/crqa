\name{CTcrqa}
\alias{CTcrqa}

\title{Contingency Table Cross-Recurrence Quantification Analysis}

\description{
  Recurrence is calculated by means of contingency tables,
  and it can only be used on categorical time-series. First,
  it finds the common state, or categories, shared by the
  two-times series, then it builds up a contingency table
  counting the co-occurrences of stateA-stateB between the
  two-series. The diagonal of the CT is where the
  recurrence profile is calculated, as along the diagonal, 
  the states are identical.
}

\usage{ CTcrqa(ts1, ts2, par) }

\arguments{
  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{par}{A list of arguments to pass to the function:
    datatype = a string specifying whether the time-series
    is 'numerical' or 'categorical'.
    thrshd = a constant indicating the maximum difference
    between time-series lenghts that is tollerated.
    lags = a numerical vector for the delays,
    e.g., seq(1,100, 1)}
}


\value{A cross-recurrence profile of the two time-series
with length equal to the number of delays considered}

\references{

  Dale, R., Warlaumont, A. S. and Richardson, D. C. (2011).
  Nominal cross recurrence as a generalized lag sequential
  analysis for behavioral streams. International Journal of
  Bifurcation and Chaos, 21, 1153-1161.
  (special issue on recurrence)

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\note{
  For every delay, a full co-occurence matrices of the states
  of the two-series is obtained. At the moment, only the
  diagonal is used. However, the function could also be used
  to track the probability of encountering states for the
  different delays.
}


\seealso{\code{\link{calcphi}}}

\examples{

## simulate two dichotomous series
tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 50)
ts1 = tS[1,]; ts2 = tS[2,]

## check data(crqa) for alternative data
## (e.g., RDts1, RDts2)

par = list(lags = seq(1, 40, 1), datatype = "categorical", thrshd = 8);
res = CTcrqa(ts1, ts2, par)

## show profile

plot(seq(1,length(res),1), res, xlab = "Delays",
ylab = "Recurrence", type = "l", lwd = 3)

}