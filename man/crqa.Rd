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
  are extracted (explained below). Since version 2.1.0 the
  function also exposes an approximative path
  (\code{method = "aRQA"}) that scales to very long series.}
\usage{
crqa(ts1, ts2, delay = 1, embed = 1, rescale = 0,
     radius = 0.001, normalize = 0, mindiagline = 2, minvertline = 2,
     tw = 0, whiteline = FALSE, recpt = FALSE, side = "both",
     method = "rqa", metric = "euclidean", datatype = "continuous",
     rr_denom = "full")
}
\arguments{
  \item{ts1}{First time-series dataset.}
  \item{ts2}{Second time-series dataset.}
  \item{delay}{The delay unit by which the series are lagged.}
  \item{embed}{The number of embedding dimension for
    phase-reconstruction, i.e., the lag intervals.}
  \item{rescale}{Rescale the distance matrix;
    if rescale = 0 (do nothing);
    if rescale = 1 (mean distance of entire matrix);
    if rescale = 2 (maximum distance of entire matrix);
    if rescale = 3 (minimum distance of entire matrix);
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
  \item{tw}{The Theiler window parameter (default 0).
    Width of the diagonal band around the line of identity
    that is blanked from the recurrence matrix.}
  \item{whiteline}{Logical; if TRUE, also compute statistics on
    the white vertical lines (gaps between recurrent stretches
    in the same column). Three additional fields are returned
    in this case: \code{wmean}, \code{wmax}, \code{wENTR}.
    Since v2.1.0 these are extracted from the same sparse-index
    pass that produces the black-line statistics, so enabling
    \code{whiteline} carries essentially no extra cost.}
  \item{recpt}{A logical flag indicating whether
    measures of cross-recurrence are calculated directly
    from a recurrent plot (TRUE) or not (FALSE).}
  \item{side}{A string indicating whether recurrence measures
    should be calculated in the 'upper' triangle of the RP,
    the 'lower' triangle, or 'both'. The line of identity is
    automatically excluded for 'upper' and 'lower'.}
  \item{method}{A string indicating the type of recurrence
    analysis to perform. Options:
    \code{"rqa"} (auto-recurrence on a single series),
    \code{"crqa"} (cross-recurrence on two series),
    \code{"mdcrqa"} (multidimensional cross-recurrence on two
    multidimensional series), or \code{"aRQA"} (since v2.1.0:
    approximative auto-RQA via phase-space histogram binning,
    Schultz et al. 2015 Phys. Lett. A 379:997-1011 and Spiegel
    et al. 2016). The aRQA path scales to N >> 10^5 (where
    \code{"rqa"}/\code{"crqa"} run out of memory) at the cost
    of a small approximation error: a few percent on DET for
    stochastic data, larger bias for strongly deterministic
    systems. Only \code{RR}, \code{DET} and \code{L} are
    populated for aRQA; \code{RP}, line-length distributions
    and vertical-line measures are returned as \code{NA}.}
  \item{metric}{A string to indicate the type of
    distance metric used (default \code{"euclidean"}).
    See \code{\link[rdist]{cdist}} for the full list of
    supported metrics.}
  \item{datatype}{A string (\code{"continuous"} or
    \code{"categorical"}) indicating the nature of the data.}
  \item{rr_denom}{Choice of denominator for the recurrence rate
    (new in v2.1.0). The literature does not settle on a single
    convention when a Theiler window or one-sided mask is in use:
    \describe{
      \item{\code{"full"} (default)}{Denominator is \code{N1 * N2}
        (all cells in the recurrence plot), regardless of how many
        were blanked by the Theiler window or side mask.
        Formula: \code{RR = 100 * numrecurs / (N1 * N2)}.
        Reproduces the historical convention of \code{crqa} <=
        2.0.7, of the original Marwan (2007, Phys. Rep. 438)
        \code{RR = (1/N^2) Sum R_ij} definition, and of the
        TOCSY MATLAB toolbox. Use this for apples-to-apples
        comparisons with other RQA tools or with published
        analyses that calibrated \code{radius} against this
        convention.}
      \item{\code{"valid"}}{Denominator counts only the cells that
        survived the Theiler band and the side mask.
        Formula: \code{RR = 100 * numrecurs / (N1 * N2 - theiler_exclusion(N1, N2, tw))}.
        \code{RR} then reflects the fraction of \emph{analysable}
        cell pairs that are recurrent, and is internally consistent
        with how \code{DET} and \code{ENTR} are already computed
        (their denominator is \code{numrecurs}, the post-Theiler
        count). The excluded cell count is given by the new
        \code{\link{theiler_exclusion}} helper, which works for
        both square and rectangular RPs.}
    }
    For \code{tw == 0} and \code{side == "both"} both
    conventions give the same value, so this argument matters
    only when one of those is non-default.}
}

\details{
  Version 2.1.0 modernised \code{crqa()} along four axes,
  preserving exact numerical identity with the previous CRAN
  release (v2.0.7) on every shared measure:

  \enumerate{
    \item Sparse-index line scan
      (\code{line_stats()}, internal). Replaces the
      \code{spdiags()}-based diagonal extraction. Same outputs,
      ~5x faster for typical N.
    \item In-place rescale. When \code{rescale > 0} the
      rescaling factor is folded into \code{radius}
      (\code{dm / s <= r} <=> \code{dm <= r * s}) so the N*N
      distance matrix is no longer copied. Halves peak RAM at
      \code{rescale > 0}.
    \item Fused C++ inner loop (when \code{metric} is
      \code{"euclidean"}, \code{"maximum"} or
      \code{"manhattan"} -- the typical case). Distance,
      threshold and line statistics are computed in a single
      pass that never materialises \code{dm} or the dense
      recurrence matrix. Memory is O(N + nnz) instead of
      O(N^2); runtime is 5-15x faster than the legacy path at
      N >= 2000.
    \item OpenMP parallelism (since v2.1.0 Stage 3c). The
      fused kernel's distance-and-threshold loop is parallelised
      across post-transpose columns. Results are bit-identical
      to the serial reference at any thread count -- only
      floating-point comparisons happen in the parallel region,
      no reductions. crqa respects the standard
      \code{OMP_NUM_THREADS} environment variable, which must
      be set in the shell BEFORE launching R (libgomp caches
      its thread count at first parallel region, so
      \code{Sys.setenv()} from inside R is too late). When
      OpenMP is unavailable (some macOS source builds with
      Apple clang lacking libomp) the kernel transparently
      falls back to single-threaded execution. Typical speedup
      at N >= 10000 with sparse recurrence: ~2x on 4 cores,
      ~3-4x on 8 cores.
  }

  The default \code{rr_denom = "full"} preserves the v2.0.7 RR
  convention for backwards compatibility. Pass
  \code{rr_denom = "valid"} for a denominator that excludes
  Theiler-blanked cells, making RR internally consistent with
  DET and ENTR. When \code{side != "both"} or \code{tw > 0},
  the two options diverge; consider explicitly setting
  \code{rr_denom} in those cases.

  Other metrics (\code{"canberra"}, \code{"minkowski"} etc.)
  still go through the legacy \code{cdist()}-based path. All
  metrics produce identical results to CRAN v2.0.7 on the
  shared measure set.
}

\value{
  If a recurrence plot (RP) can be calculated and recurrence
  is observed, the function returns a list of measures.
  Otherwise the values are 0 or \code{NA}.
  \item{RR}{Percentage of recurrent points (range 0-100). The
    denominator depends on \code{rr_denom} (see above).}
  \item{DET}{Proportion of recurrent points forming diagonal
    line structures of length >= \code{mindiagline}.}
  \item{NRLINE}{Total number of diagonal lines >= \code{mindiagline}.}
  \item{maxL}{Length of the longest diagonal line segment
    (excluding the main diagonal if \code{tw > 0}).}
  \item{L}{Average length of diagonal lines >= \code{mindiagline}.}
  \item{ENTR}{Shannon entropy of the distribution of diagonal
    line lengths >= \code{mindiagline}.}
  \item{rENTR}{\code{ENTR} normalised by the number of distinct
    line lengths observed.}
  \item{LAM}{Proportion of recurrent points forming vertical
    line structures of length >= \code{minvertline}.}
  \item{TT}{Trapping Time: average length of vertical lines
    >= \code{minvertline}.}
  \item{catH}{Entropy of categorical recurrence plots based
    on rectangular block structures (only when
    \code{datatype = "categorical"}, \code{side = "both"} and
    \code{radius <= 0.1}).}
  \item{max_vertlength}{Maximum vertical line length.}
  \item{wmean, wmax, wENTR}{White vertical line statistics
    (mean / max length, and Shannon entropy of their length
    distribution = recurrence-time entropy in the sense of
    Faure & Korn 2004 Phys. Lett. A 332:329-339 and Little et
    al. 2007 Biomed. Eng. Online 6:23). Returned as \code{NA}
    when \code{whiteline = FALSE}. New in v2.1.0.}
  \item{RP}{The recurrence plot as a sparse matrix
    (\code{\link[Matrix]{Matrix}}). \code{NA} when
    \code{method = "aRQA"} (the matrix is never built on that
    path; see \code{\link{plot_rp}} for the user-facing error).}
}

\references{
Coco, M. I., and Dale, R. (2014). Cross-recurrence quantification analysis of categorical and continuous time series: an R package. Frontiers in psychology, 5, 510.

Wallot, S. (2018). Multidimensional Cross-Recurrence Quantification Analysis (MdCRQA) - a method for quantifying correlation between multivariate time-series. Multivariate behavioral research.

Schultz, D., Spiegel, S., Marwan, N., Albayrak, S. (2015). Approximation of diagonal line based measures in recurrence quantification analysis. Phys. Lett. A 379(14-15):997-1011.

Marwan, N., Romano, M.C., Thiel, M., Kurths, J. (2007). Recurrence plots for the analysis of complex systems. Phys. Rep. 438(5-6):237-329.}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\note{Original bits of this code were translated from a Matlab
  version provided by Rick Dale, and created during
  the Non-Linear Methods for Psychological Science
  summer school held at the University of Cincinnati
  in 2012. The multi-dimensional method for \code{crqa} was
  developed together with Sebastian Wallot.}

\seealso{
\code{\link{theiler_exclusion}},
\code{\link{aRQA}}, \code{\link{spdiags}}, \code{\link{simts}}
}

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

# print everything except the recurrence plot itself
print(ans[!names(ans) \%in\% "RP"])

}

\keyword{ts}
