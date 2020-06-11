\name{windowdrp}
\alias{windowdrp}

\title{Windowed Recurrence Profile}

\description{
  
  A recurrence plot is computed in overlapping
  windows of a specified size for a number of
  delays smaller than the size of the window.
  In every window, the recurrence value for the
  different delays is calculated.
  A mean is then taken across the delays to obtain
  a recurrence value in that particular window.

}

\usage{
windowdrp(ts1, ts2, windowstep, windowsize, lagwidth,
radius, delay, embed, rescale, normalize, mindiagline, minvertline,
tw, whiteline, recpt, side, method, metric, datatype)

}

\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{windowstep}{Interval by which the window is
    moved.}
  \item{windowsize}{The size of the window}  
  \item{lagwidth}{The number of delays to be
    considered within the window}
  \item{radius}{For numeric time-series, the
    cutoff distance to accept or reject
    two-points as recurrent}
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
  
  It returns a list of arguments where:
  \item{profile}{Time-course windowed recurrence profile}
  \item{maxrec}{Maximal recurrence observed along the time-course}
  \item{maxlag}{The point where maximal recurrence is observed}

}

\references{

  Boker, S. M., Rotondo, J. L., Xu, M., and King, K. (2002).
  Windowed cross-correlation and peak picking
  for the analysis of variability in the
  association between behavioral time series.
  Psychological Methods, 7(3), 338.

}
\author{Moreno I. Coco (moreno.cocoi@gmail.com)
  and Rick Dale (rdale@ucmerced.edu)}


\seealso{\code{\link{drpfromts}}}

\examples{

# use the available data
data(crqa) 

listener = eyemovement$listener
narrator = eyemovement$narrator

delay = 1; embed = 1; rescale = 1; radius = 0.001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = 'crqa'; metric = 'euclidean';  
datatype = "continuous"; windowsize =  200; 
lagwidth = 40; windowstep = 10

ans = windowdrp(narrator, listener, windowstep, windowsize, lagwidth, 
                radius, delay, embed, rescale, normalize, 
                mindiagline, minvertline, tw, 
                whiteline, recpt, side, method, metric, 
                datatype)


profile = ans$profile; maxrec = ans$maxrec; maxlag = ans$maxlag

plot(profile, type = 'l')

}