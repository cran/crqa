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

wincrqa(x, y, step, windowsize, lagwidth, delay, embed, rescale,
radius, normalize, mindiagline, minvertline, tw, whiteline, recpt)

}


\arguments{

  \item{x}{First time-series}
  \item{y}{Second time-series}
  \item{step}{Interval by which the window is moved.}
  \item{windowsize}{The size of the window}  
  \item{lagwidth}{The number of delays to be considered}
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
  \item{recpt}{A logical flag indicating whether
    measures of cross-recurrence are calculated directly
    from a recurrent plot (TRUE) or not (FALSE).}
  
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
  about it will be warned.}


\seealso{\code{\link{crqa}}}

\examples{

## simulate two dichotomous series
tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 200)
ts1 = tS[1,]; ts2 = tS[2,]

## check data(crqa) for alternative data
## (e.g., RDts1, RDts2)

step =  10; windowsize =  100; lagwidth = 20;
delay = 1; embed =  1 ; rescale =  1;
radius = 0.00001; normalize = 0;
minvertline = 2; mindiagline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE

## it returns matrix with all measures for the
## different windows where values are found

res = wincrqa(ts1, ts2, step, windowsize,
lagwidth, delay, embed, rescale, radius, normalize,
mindiagline, minvertline, tw, whiteline, recpt)

str(res)

}
