\name{crqa}
\alias{crqa}
\title{Cross recurrence measures of two time-series,
  time-delayed and embedded in higher dimensional space
}
\description{Core cross recurrence function, which examines
  recurrent structures between time-series, which are
  time-delayed and embedded in higher dimensional space.
  The approach compares the phase space trajectories of two
  time-series in the same phase-space when delays are
  introduced. A distance matrix between the two-series,
  delayed and embedded is calculated. Several measures
  representative of the interaction between two series are
  extracted (explained below).}
\usage{
crqa(ts1, ts2, delay, embed, rescale, radius, normalize,
mindiagline, minvertline, tw, whiteline, recpt, side, checkl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ts1}{First time-series.}
  \item{ts2}{Second time-series.}
  \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
    if rescale = 0 (do nothing);
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
  \item{checkl}{A list with four arguments:
      do = TRUE|FALSE; normalize (or not) the length of ts
      if do == TRUE, then the arguments of checkts() needs to be passed.
      datatype = (numerical, categorical) - nature of ts
      thrshd: number of timepoints we tollerate difference between ts
      pad: the two series have to be padded (TRUE) or chopped (FALSE)}	
}

\details{
  We recommend setting whiteline = FALSE, as the current
  version of the library does not make use of such
  information to extract measures of cross-recurrence.
}

\value{
  If CRQA can be calculated and recurrence is found, it returns
  a list with different measures extracted from the recurrence
  plot. Otherwise, the values for the output arguments will be
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
  \item{RP}{The Recurrence Plot sparse matrix data}
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\note{Part of this code was translated from a Matlab
  version provided by Rick Dale, and created during 
  the Non-Linear Methods for Psychological Science
  summer school held at the University of Cincinnati
  in 2012}

\seealso{ \code{\link{tt}}, \code{\link{checkts}},
  \code{\link{spdiags}}, \code{\link{simts}}}

\examples{

## simulate two dichotomous series
tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 100)
ts1 = tS[1,]; ts2 = tS[2,]

## check data(crqa) for alternative data
## (e.g., RDts1, RDts2)

## initialize the parameters

delay = 1; embed = 1; rescale = 1; radius = 0.001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
checkl = list(do = FALSE, thrshd = 3, datatype = "categorical",
    pad = TRUE)

ans = crqa(ts2, ts1, delay, embed, rescale, radius, normalize,
mindiagline, minvertline, tw, whiteline, recpt, side, checkl)

print(ans[1:9]) ## last argument of list is the cross-recurrence plot
RP = ans$RP ## take out RP

}

\keyword{ts}
