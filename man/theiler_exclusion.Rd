\name{theiler_exclusion}
\alias{theiler_exclusion}
\title{Number of cells excluded by a Theiler window}
\description{
Counts the cells of an \code{m x n} recurrence-plot matrix that fall
inside a Theiler band of half-width \code{w} around the line of
identity. Useful when computing a recurrence rate that excludes the
Theiler-blanked cells from the denominator (as the new
\code{rr_denom = "valid"} mode of \code{\link{crqa}} does).
}
\usage{
theiler_exclusion(m, n = m, w = 1L)
}
\arguments{
  \item{m}{Number of rows of the recurrence matrix.}
  \item{n}{Number of columns. Defaults to \code{m} (square matrix).}
  \item{w}{Theiler-window half-width:
    \code{w = 0} means no exclusion (returns 0);
    \code{w = 1} excludes the main diagonal only (returns
    \code{min(m, n)});
    \code{w = k} excludes the \code{2k - 1} diagonals centred on the
    line of identity.}
}
\details{
For square matrices the count reduces to the closed form
\code{m * (2 * w - 1) - w * (w - 1)} given by pjbruna's community PR
on the upstream \pkg{crqa} repository. For rectangular matrices the
asymmetric diagonal lengths are summed exactly, so the function is
correct in either case.
}
\value{Integer count of excluded cells.}

\seealso{\code{\link{crqa}}}

\examples{
theiler_exclusion(100, 100, 1)   # main diagonal:           100
theiler_exclusion(100, 100, 2)   # main +/- 1 diagonal:     298
theiler_exclusion(3, 5, 2)       # rectangular case:          8
theiler_exclusion(100, 100, 0)   # no exclusion:              0
}

\keyword{utilities}
