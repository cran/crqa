\name{piecewiseRQA}
\alias{piecewiseRQA}
\title{Compute recurrence plots for long time-series data series
using a block (piece-wise) method.}

\description{This is a convenience function which breaks down 
the computation of large recurrence plots into a collection of
smaller recurrence plots. It can ease speed and memory issues
if an appropriate size for the block is found.}
  
\usage{
piecewiseRQA(ts1, ts2, blockSize, delay = 1, embed = 1, rescale = 0,
radius = 0.001, normalize = 0, mindiagline = 2, minvertline = 2,
tw = 0, whiteline = FALSE, recpt = FALSE, side = "both",
method = "crqa", metric = "euclidean", datatype = "continuous",
typeRQA = "full", windowsize = NA,
workers = 1L, rr_denom = "full")
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
    if normalize =  2 (z-score).
    Since v2.1.0, normalization is applied once to the \emph{full} series
    before block decomposition, so all blocks share the same scale and
    the assembled recurrence plot is internally consistent.}
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
\item{workers}{Integer; number of parallel workers (default 1 =
  serial). New in v2.1.0; see \code{\link{wincrqa}} for the
  parallel-execution caveats (including the OpenMP interaction).}
\item{rr_denom}{Forwarded to \code{\link{crqa}}; choice of
  recurrence-rate denominator. New in v2.1.0; see
  \code{\link{crqa}} for details.}

}

\details{
  \code{piecewiseRQA} is a structural tool: it divides the time series
  into consecutive non-overlapping blocks and computes a separate
  recurrence plot for each block, then aggregates the measures. The
  primary use case is decomposing a long recording into epochs and
  examining how recurrence structure changes across them.

  Note: since v2.1.0 the core \code{crqa()} function uses a fused C++
  inner loop that never materialises the full N x N distance matrix
  (memory is O(N + nnz)), so \code{piecewiseRQA} is no longer needed
  as a memory or speed workaround. For analysing a single long series
  consider \code{crqa()} directly (exact measures) or
  \code{crqa(\ldots, method = "aRQA")} (O(N) memory approximation).

  Set \code{workers > 1} to process blocks in parallel via the
  \pkg{future} + \pkg{furrr} backends.
}

\value{
  If an RP can be calculated and recurrence is found, 
  the piecewiseRQA will return exactly the same measures
  as crqa() if the typeRQA is set to 'full' and
  drpfromts() if the typeRQA is set to 'diagonal'.
  Please refer to the help file for those two functions 
  for details about the measures.
  \item{RP}{The Recurrence Plot sparse matrix data}

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com) based 
on Matlab code by Sebastian Wallot}

\seealso{\code{\link{crqa}}, \code{\link{spdiags}}, \code{\link{simts}}}

\examples{

\donttest{

data(crqa)

listener = eyemovement$listener
narrator = eyemovement$narrator

blockSize = 200; delay = 1; embed = 1; rescale = 0; radius = 0.001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = "crqa"; metric = "euclidean"; datatype = "categorical"
typeRQA = "full"; windowsize = NA

pieceRP = piecewiseRQA(narrator, listener, blockSize, delay, embed, rescale,
                       radius, normalize, mindiagline, minvertline,
                       tw, whiteline, recpt, side,
                       method, metric, datatype, typeRQA,
                       windowsize)

## print recurrence measures across blocks (excluding the RP objects)
print(pieceRP[!names(pieceRP) \%in\% "RP"])

}

}

\keyword{ts}
