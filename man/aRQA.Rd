\name{aRQA}
\alias{aRQA}
\title{Approximative Recurrence Quantification Analysis}
\description{
Computes approximate recurrence rate (RR), diagonal-line density (DET)
and mean diagonal-line length (L) for one or two (auto- or cross-)
embedded time-series via phase-space histogram binning, without ever
materialising the recurrence matrix. Run-time is essentially linear in
the number of embedded points; the algorithm is the one published in
Schultz et al. (2015) and Spiegel et al. (2016) and is the practical
route for \code{N >> 10^5}, where the exact \code{\link{crqa}} pipeline
runs out of memory.

\strong{Note:} The returned \code{DET} is an approximated measure of
diagonal-line density derived from histogram statistics, not the exact
determinism obtained by scanning the recurrence matrix. Treat it
accordingly when interpreting results, especially for strongly
deterministic dynamics (see Details).

This function is also reachable as \code{crqa(..., method = "aRQA")}.
}
\usage{
aRQA(ts1, ts2 = ts1, delay = 1, embed = 1, radius = 0.1,
     mindiagline = 2, normalize = 0)
}
\arguments{
  \item{ts1, ts2}{Numeric vectors (auto-RQA: \code{ts2 = ts1};
    cross-RQA: distinct vectors). For multidimensional data pass
    matrices with rows = time, cols = dimensions.}
  \item{delay}{Embedding delay tau (>= 1).}
  \item{embed}{Embedding dimension m (>= 1).}
  \item{radius}{Threshold epsilon. Used as the side length of the
    hypercubes that partition phase space.}
  \item{mindiagline}{Minimum diagonal line length (used as L_min in
    Schultz et al.'s closed-form expressions for DET and L).}
  \item{normalize}{0 = none (default), 1 = unit-interval rescale,
    2 = z-score. Same convention as \code{\link{crqa}}.}
}

\details{
The algorithm implements Eqs. (15-16) of Schultz, Spiegel, Marwan &
Albayrak (2015): partition the m-dimensional embedded phase space into
hypercubes of side \code{radius}, build the histogram h of cube
occupancies, and approximate the recurrence rate as
\deqn{RR \approx \frac{1}{N^2}\sum_b h(b)^2.}

For diagonal-line measures, RR is computed at three "stack heights"
(L = 1, mindiagline, mindiagline + 1) by histogramming L-grams of
consecutive embedded points and using the closed-form relations
\deqn{DET(L_{min}) \approx \frac{L_{min} \cdot RR(L_{min}) - (L_{min}-1) \cdot RR(L_{min}+1)}{RR(1)},}
\deqn{L(L_{min})   \approx \frac{L_{min} \cdot RR(L_{min}) - (L_{min}-1) \cdot RR(L_{min}+1)}{RR(L_{min}) - RR(L_{min}+1)}.}

\strong{Approximation error.} Two embedded points within \code{radius}
in max-norm distance can fall in adjacent hypercubes (false negatives
across bin boundaries). Errors are typically a few percent on \code{DET}
for stochastic / weakly deterministic data (the regime targeted by the
original aRQA papers) and grow for strongly deterministic dynamics
where most recurrent points lie on long diagonals (the L > 1 histograms
become sparse — the "curse of dimensionality" hits the L-stacks). For
publication-quality RQA on strongly deterministic systems, fall back to
\code{method = "rqa"} or \code{"crqa"}.

\strong{Short series.} aRQA's advantage is memory and speed at large N.
For series shorter than roughly 1,000 embedded points the exact
\code{\link{crqa}} pipeline (\code{method = "rqa"} or \code{"crqa"})
finishes in milliseconds and returns precise values. Using
\code{method = "aRQA"} on short data gains nothing and introduces
unnecessary approximation error; a warning is issued in that case.

\strong{Parameters silently ignored.} Because the algorithm never
materialises the recurrence matrix, the following \code{\link{crqa}}
parameters have no effect on the aRQA path: \code{tw}, \code{side},
\code{whiteline}, \code{minvertline}, \code{metric}, \code{recpt},
\code{rescale}.
}

\value{
A list matching the structure of \code{\link{crqa}}'s return value, with
\code{NA} for measures that cannot be derived from the histogram alone:
\item{RR}{Approximate recurrence rate (percentage).}
\item{DET}{Approximate diagonal-line density (percentage): the estimated
  proportion of recurrent points that belong to diagonal structures of
  length >= \code{mindiagline}. Derived from histogram statistics, not
  from scanning the recurrence matrix.}
\item{L}{Approximate mean diagonal-line length.}
\item{NRLINE, maxL, ENTR, rENTR, LAM, TT, max_vertlength, catH}{\code{NA}}
\item{RP}{\code{NA} — the recurrence matrix is never built}
}

\references{
Schultz, D., Spiegel, S., Marwan, N., Albayrak, S. (2015). Approximation of diagonal line based measures in recurrence quantification analysis. Phys. Lett. A 379(14-15):997-1011.

Spiegel, S., Schultz, D., Marwan, N. (2016). Approximate Recurrence Quantification Analysis (aRQA) in Code of Best Practice. In: Webber, Ioana, Marwan (eds), Recurrence Plots and Their Quantifications: Expanding Horizons. Springer.
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com); algorithm by
Schultz, Spiegel, Marwan & Albayrak (2015).}

\seealso{\code{\link{crqa}}}

\examples{
set.seed(1)
x <- as.numeric(arima.sim(model = list(ar = 0.7), n = 2000))
res <- aRQA(x, x, delay = 1, embed = 3, radius = 0.5,
            mindiagline = 2, normalize = 2)
res$RR
res$DET
res$L
}

\keyword{ts}
