\name{crqa}
\alias{crqa}
\title{Auto, cross and multidimensional recurrence measures 
of one, two or multiple time-series, time-delayed and embedded 
in higher dimensional space
}
\description{Core recurrence function, which examines
  recurrent structures of a single (rqa), two (crqa), 
  or multidimensional time-series (mdcrqa), which are
  time-delayed and embedded in higher dimensional space.
  The approach compares the phase space trajectories of 
  the time-series in the same phase-space when delays are
  introduced. A distance matrix between the time-series,
  delayed and embedded is calculated. Several measures
  representative of the underlying dynamics of the system
  are extracted (explained below).}
\usage{
crqa(ts1, ts2, delay, embed, rescale, radius, normalize,
mindiagline, minvertline, tw, whiteline, recpt, side, method,
metric, datatype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ts1}{First time-series dataset.}
  \item{ts2}{Second time-series dataset}
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
}

\details{
  We recommend setting whiteline = FALSE, as the current
  version of the library does not make use of such
  information to extract recurrence measures.
}

\value{
  If a recurrence plot (RP) can be calculated and hence 
  recurrence observed the function will returna a list with 
  different measures extracted. 
  Otherwise, the values for the output arguments will be
  either 0 or NA.
  \item{RR}{The percentage of recurrent points falling within
    the specified radius (range between 0 and 100)}
  \item{DET}{Proportion of recurrent points forming diagonal
    line structures.}
  \item{NRLINE}{The total number of lines in the recurrent plot}
  \item{maxL}{The length of the longest diagonal line
    segment in the plot, excluding the main diagonal}
  \item{L}{The average length of line structures}
  \item{ENTR}{Shannon information entropy of
    diagonal line lengths longer than the minimum length}
  \item{rENTR}{Entropy measure normalized by the number of
    lines observed in the plot. Handy to compare across contexts
    and conditions}
  \item{LAM}{Proportion of recurrent points forming vertical
    line structures}
  \item{TT}{The average length of vertical line structures}
  \item{catH}{Entropy of categorical recurrence plots based on rectangular block structures}
  \item{RP}{The Recurrence Plot sparse matrix data}
}

\references{
Coco, M. I., and Dale, R. (2014). Cross-recurrence quantification analysis of categorical and continuous time series: an R package. Frontiers in psychology, 5, 510.

Wallot, S. (2018). Multidimensional Cross-Recurrence Quantification Analysis (MdCRQA) a method for quantifying correlation between multivariate time-series. Multivariate behavioral research, 1-19}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\note{Original bits of this code were translated from a Matlab
  version provided by Rick Dale, and created during 
  the Non-Linear Methods for Psychological Science
  summer school held at the University of Cincinnati
  in 2012. The multi-dimensional method for the crqa function 
  has been written together with Sebastian Wallot 
  (sebastian.wallot at aesthetics.mpg.de )}

\seealso{\code{\link{spdiags}}, \code{\link{simts}}}

\examples{

# use the available data
data(crqa) 

listener = eyemovement$listener
narrator = eyemovement$narrator

delay = 1; embed = 1; rescale = 0; radius = .1;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = 'crqa'; metric = 'euclidean';  
datatype = "categorical"

ans = crqa(narrator, listener, delay, embed, rescale, radius, normalize, 
           mindiagline, minvertline, tw, whiteline, recpt, side, method,
           metric, datatype)

print(ans[1:10]) ## last argument of list is the recurrence plot

}

\keyword{ts}
