\name{crqa-package}
\alias{crqa-package}

\docType{package}

\title{
  Cross-Recurrence Quantification Analysis for Continuous and
  Categorial Time-series
}

\description{
  CRQA is a package to perform cross-recurrence quantification
  analysis between two time-series, of either categorical or
  continuous values. It provides different methods
  for profiling cross-recurrence, i.e., only looking
  at the diagonal recurrent points, as well as more
  in-depth measures of the whole cross-recurrence plot,
  e.g., percentage recurrence.
}
\details{
\tabular{ll}{
Package: \tab crqa\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-20\cr
License: \tab GPL >= 2\cr
}


\code{crqa}: core cross recurrence function, which examines
recurrent structures between time-series, which are
time-delayed and embedded in higher dimensional space.

\code{drpdfromts}: cross-recurrence profile of two-times series
  
\code{CTcrqa}: calculates cross-recurrence of categorical
time-series by means of Contingency Tables.
Usefull if co-occurences between states has to be tracked.

\code{optimizeParam}: find the optimal radius, delay,
and embedding dimension parameters to compute cross-recurrence
of two continuous time-series.

\code{runcrqa}: Wrapper calling different methods to
compute CRQA. For each method, an appropriate set of
parameters should be chosen.

\code{windowdrp}: windowed-recurrence profile, similar to window
cross-correlation analysis, it calculates how recurrence
of the two time-series develop over time. 

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}

\references{
  
  Webber Jr, C. L., and Zbilut, J. P. (2005). Recurrence
  quantification analysis of nonlinear dynamical systems.
  Tutorials in contemporary nonlinear methods for the
  behavioral sciences, 26-94.
  Marwan, N., and Kurths, J. Nonlinear analysis of bivariate
  data with cross recurrence plots. Physics Letters A 302.5
  (2002): 299-307.
  
}

\keyword{ package }

\examples{

data(crqa) ## load data 

## cross-recurrence diagonal profile
drpdfromts(RDts1, RDts2, ws = 40, datatype = "categorical")


}
