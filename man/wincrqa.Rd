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
radius = 0.001, rescale = 0, normalize = 0,
mindiagline = 2, minvertline = 2,
tw = 0, whiteline = FALSE, recpt = FALSE, side = "both",
method = "crqa", metric = "euclidean", datatype = "continuous",
trend = FALSE, workers = 1L,
rr_denom = "full")

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
    if normalize =  2 (z-score).
    Since v2.1.0, normalization is applied once to the \emph{full} series
    before windowing, so all windows share the same scale.
    In earlier versions it was applied independently per window.}
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
\item{workers}{Integer; number of parallel workers (default 1 =
  serial). Set greater than 1 to dispatch the per-window
  \code{crqa()} calls via \pkg{future} + \pkg{furrr}. Note: spawning
  workers has a fixed startup cost of a few seconds; parallel
  execution only pays off when individual windows are themselves
  expensive (large \code{N} per window or many windows). For typical
  exploratory use (short series, modest number of windows) the
  serial default is faster. New in v2.1.0.

  \strong{Interaction with OpenMP.} Since v2.1.0 the fused C++
  kernel inside \code{crqa()} is itself OpenMP-parallel and by
  default uses all available cores. Combining \code{workers > 1}
  with the default OpenMP setting can over-subscribe the CPU
  (e.g. \code{workers = 4} on a 4-core machine would spawn up to
  16 threads). When using \code{workers > 1}, set
  \code{OMP_NUM_THREADS=1} in the shell environment BEFORE
  launching R, so each worker uses a single OpenMP thread and
  the parallelism comes entirely from \code{workers}.
  \code{Sys.setenv()} from inside R is too late: libgomp caches
  the thread count at first parallel region.}
\item{rr_denom}{Forwarded to \code{\link{crqa}}; choice of
  recurrence-rate denominator. \code{"full"} (default) uses all
  \code{N1 * N2} cells and reproduces the historical convention;
  \code{"valid"} excludes Theiler/side-blanked cells from the
  denominator. New in v2.1.0; see \code{\link{crqa}} for details.}

}


\value{

  It returns a matrix where the rows are the
  different windows explored, and the columns
  are the recurrence measures observed
  in that particular window.
  Refer to \code{crqa} for the values returned.


}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)
Alexandra Paxton \email{alexandra.paxton@uconn.edu}}

\note{If no-recurrence is found in a window,
  that window will not be saved, and a message
  about it will be warned.
  TREND is implemented following a solution proposed
  by Norbert Marwan, and translated here in R,
  for those who have asked him.
  He, however, warns that this measure might strongly depend
  on the chosen settings to calculate crq.
  Relying blindly on such measure can, therefore, produce misleading results.
  Also, we enabled the possibility to input directly a RP with recpt = T, 
  and so extract the measures. This implies that it will not be the same 
  as by conducting phase-space reconstruction on subsets of the time series.  
  So, please, make sure why you are doing it and why.
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
windowsize = 200; windowstep = 100
trend = FALSE

ans = wincrqa(listener, narrator, windowstep, windowsize, delay, embed,
                    radius, rescale, normalize, mindiagline, minvertline,
                    tw, whiteline, recpt, side, method, metric, 
                    datatype, trend)

## other recurrence measures are available in ans
profile = as.numeric(ans$RR) 

plot(profile, type = 'l')

}
