\name{wincrqa}
\alias{wincrqa}

\title{Window Cross-Recurrence Measures}

\description{

  It computes cross-recurrence is calculated in overlapping
  windows of a certain size for a number of delays smaller
  than the size of the window. For each window,
  a cross-recurrence plot is build and measures of it extracted.
  
}

\usage{

wincrqa(ts1, ts2, windowstep, windowsize, delay, embed, rescale,
radius, normalize, mindiagline, minvertline, tw, whiteline, trend)

}


\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{windowstep}{Interval by which the window is moved.}
  \item{windowsize}{The size of the window}  
  \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
    if rescale = 1 (mean distance of entire matrix);
    if rescale = 2 (maximum distance of entire matrix).}
  \item{radius}{A threshold, cut-off, constant used to
      decide whether two points are recurrent or not.}
  \item{normalize}{Normalize the time-series;
    if normalize = 0 (do nothing);
    if normalize = 1 (Unit interval);
    if normalize =  2 (z-score).}
  \item{mindiagline}{A minimum diagonal length of
    recurrent points. Usually set to 2, as it takes
    a minimum of two points to define any line.}
  \item{minvertline}{A minimum vertical length of
    recurrent points.}
  \item{tw}{The size of the Theiler window}
  \item{whiteline}{A logical flag to calculate (TRUE)
    or not (FALSE) empty vertical lines.}
  \item{trend}{A logical flag indicating whether
    the TREND should be computed}
  
}


\value{

  It returns a matrix where the rows are the
  different windows explored, and the columns
  are the cross-recurrence measures observed
  in that particular window.
  Refer to \code{crqa} for the values returned.


}

\author{Moreno I. Coco (moreno.cocoi@gmail.com) }

\note{If no-recurrence is found in a window,
  that window will not be saved, and a message
  about it will be warned.
  TREND is implemented following a solution proposed
  by Norbert Marwan, and translated here in R,
  for those who have asked him.
  He, however warns that this measure might strongly depend
  on the chosen settings to calculate crq.
  Relying on such measure can, therefore, produce misleading results.
  }


\seealso{\code{\link{crqa}}}

\examples{

## simulate two dichotomous series
tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 100)
ts1 = tS[1,]; ts2 = tS[2,]

## check data(crqa) for alternative data
## (e.g., RDts1, RDts2)

windowstep =  10; windowsize =  50; 
delay = 1; embed =  1 ; rescale =  1;
radius = 0.00001; normalize = 0;
minvertline = 2; mindiagline = 2;
tw = 0; whiteline = FALSE; trend = TRUE;

## it returns a list with:
## [[1]] the measures for the different windows where  values are found
## [[2]] the trend over time. 

res = wincrqa(ts1, ts2, windowstep, windowsize,
delay, embed, rescale, radius, normalize, mindiagline,
minvertline, tw, whiteline, trend)

str(res)

}
