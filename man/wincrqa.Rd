\name{wincrqa}
\alias{wincrqa}

\title{Windowed Recurrence Measures}

\description{

  A recurrence plot is computed in overlapping
  windows of a certain size for a number of delays smaller
  than the size of the window; and measures of it extracted.
  
}

\usage{

wincrqa(ts1, ts2, windowstep, windowsize, delay, embed, 
radius, rescale, normalize, mindiagline, minvertline, tw, whiteline, 
recpt, side, method, metric, datatype, trend)

}


\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{windowstep}{Interval by which the window is moved.}
  \item{windowsize}{The size of the window}  
  \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{radius}{A threshold, cut-off, constant used to
      decide whether two points are recurrent or not.}
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
\item{trend}{a boolean (TRUE or FALSE) 
  to indicate whether the TREND should be computed of the
  system}

}


\value{

  It returns a matrix where the rows are the
  different windows explored, and the columns
  are the recurrence measures observed
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

data(crqa) 

listener = eyemovement$listener
narrator = eyemovement$narrator

delay = 1; embed = 1; rescale = 0; radius = 0.001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = 'crqa'; metric = 'euclidean';  
datatype = "continuous"; 
windowsize =  100; windowstep = 20
trend = FALSE

ans = wincrqa(listener, narrator, windowstep, windowsize, delay, embed,
                    radius, rescale, normalize, mindiagline, minvertline,
                    tw, whiteline, recpt, side, method, metric, 
                    datatype, trend)

## other recurrence measures are available in ans
winres  = as.data.frame(ans$crqwin)
profile = as.numeric(winres$RR) 

plot(profile, type = 'l')

}
