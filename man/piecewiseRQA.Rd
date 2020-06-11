\name{piecewiseRQA}
\alias{piecewiseRQA}
\title{Compute recurrence plots for long time-series data series
using a block (piece-wise) method.}

\description{This is a convenience function which breaks down 
the computation of large recurrence plots into a collection of
smaller recurrence plots. It can ease speed and memory issues
if an appropriate size for the block is found.}
  
\usage{
piecewiseRQA(ts1, ts2, blockSize, delay, embed, rescale, radius,
normalize, mindiagline, minvertline, tw, whiteline, recpt, side, 
method, metric, datatype, typeRQA, windowsize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ts1}{First time-series.}
  \item{ts2}{Second time-series.}
  \item{blockSize}{The dimension of the time-series subunit in which the-recurrence plot will be computed}
   \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
    if rescale = 0 (do nothing);
    if rescale = 1 (mean distance of entire matrix);
    if rescale = 2 (maximum distance of entire matrix).
    if rescale = 3 (minimum distance of entire matrix).
    if rescale = 4 (euclidean distance of entire matrix).}
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
\item{typeRQA}{a string (full or diagonal) to indicate
whether piecewise recurrence quantification measures 
should be returned for full plot or for the diagonal profile}  
\item{windowsize}{the size of the window around the
diagonal of the recurrence (if typeRQA = diagonal)}

}

\details{
  It is important to estimate the size of the block that 
  may maximize the speed of the computation. 
  We suggest to monitor how speed and memory usage
  changes as a result of increasing the time-series and
  the block size. 
  We also recommend setting whiteline = FALSE, as the current
  version of the library does not make use of such
  information to extract measures of recurrence.
}

\value{
  If an RP can be calculated and recurrence is found, 
  the piecewiseRQA will return exactly the same measures 
  as crqa() if the typeRQA is set to 'full' and 
  drpdfromts() if the typeRQA is set to 'diagonal'.
  Please refer to the help file for those two functions 
  for details about the measures.
  \item{RP}{The Recurrence Plot sparse matrix data}

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com) based 
on Matlab code by Sebastian Wallot}

\seealso{\code{\link{crqa}}, \code{\link{spdiags}}, \code{\link{simts}}}

\examples{

## Uncomment and run locally

## generate some data using pracma

# ts1 = seq(0.1, 200, .1)
# ts1 = sin(ts1) + linspace(0, 1,length(ts1));

# ts2 = ts1

## initialize the parameters
# blockSize = 100; delay  = 15; embed  = 2; rescale = 0; radius = 0.5;
# normalize = 0; mindiagline = 2; minvertline = 2;
# tw = 1; whiteline = FALSE; recpt = FALSE; side = 'both'
# method = "crqa"; metric = 'euclidean'; datatype = "continuous"
# typeRQA = "full"; windowsize = NA

# pieceRP = piecewiseRQA(ts1, ts2, blockSize, delay, embed, rescale,
#                       radius, normalize, mindiagline, minvertline,
#                       tw, whiteline, recpt, side,
#                       method, metric, datatype, typeRQA, 
#                       windowsize)
                       
# print(unlist(pieceRP[1:10]))
                       
}

\keyword{ts}
