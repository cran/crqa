\name{crqa-package}
\alias{crqa-package}

\docType{package}

\title{
  Unidimensional and Multidimensional Methods for Recurrence Quantification Analysis
}

\description{
 Auto, Cross and Multi-dimensional recurrence quantification analysis. 
    Different methods for computing recurrence, cross vs. multidimensional
    or profile iti.e., only looking at the diagonal recurrent points, 
    as well as functions for optimization and plotting are proposed.
    in-depth measures of the whole cross-recurrence plot,
    Please refer to by Coco and Dale (2014) <doi:10.3389/fpsyg.2014.00510> 
    and Wallot (2018) <doi: 10.1080/00273171.2018.1512846>
    for further details about the method.
}

\details{

\code{crqa}: Core recurrence function, which examines
  recurrent structures of a single rqa, two crqa, 
  or multidimensional time-series mdcrqa, which are
  time-delayed and embedded in higher dimensional space.
  The approach compares the phase space trajectories of 
  the time-series in the same phase-space when delays are
  introduced. A distance matrix between the time-series,
  delayed and embedded is calculated. Several measures
  representative of the underlying dynamics of the system
  are extracted.

\code{drpfromts}: Method to explore the diagonal profile of 
  the recurrence plot (Auto, Cross, or Multi-dimensional).
  It returns the recurrence for different delays, 
  the maximal recurrence observed and the delay at which it occurred.
  
\code{lorenzattractor}: An implementation of the Lorenz dynamical system,
  which describes the motion of a possible particle, which will
  neither converge to a steady state, nor diverge to infinity;
  but rather stay in a bounded but 'chaotically' defined
  region, i.e., an attractor. 

\code{mdDelay}:Estimates time delay for embedding of a 
multi-dimensional dataset.

\code{mdFnn}: Computes the percentage of false nearest
neighbors for multidimensional time series as a function 
of embedding dimension.

\code{optimizeParam}:   Iterative procedure to examine the values
  of delay, embedding dimension and radius 
  to compute recurrence plots of one, two, 
  or more time-series.

\code{piecewiseRQA}: This is a convenience function which breaks down 
the computation of large recurrence plots into a collection of
smaller recurrence plots. It can ease speed and memory issues
if an appropriate size for the block is found.

\code{plotRP}: A convenience function to plot the RP matrix returned by 
  the crqa.

\code{simts}: A simple algorithm for producing a time-series that drives
  a second time-series using parameters, which change independent and conditional
  probability of an event to occur.

\code{wincrqa}: A recurrence plot is computed in overlapping
  windows of a certain size for a number of delays smaller
  than the size of the window; and measures of it extracted.

\code{windowdrp}:  A recurrence plot is computed in overlapping
  windows of a specified size for a number of
  delays smaller than the size of the window.
  In every window, the recurrence value for the
  different delays is calculated.
  A mean is then taken across the delays to obtain
  a recurrence value in that particular window.

}

\author{Moreno I. Coco \email{moreno.cocoi@gmail.com}
    Dan Monster \email{danm@econ.au.dk}
    Giuseppe Leonardi \email{g.leonardi@vizja.pl}
    Rick Dale \email{rdale@ucla.edu} 
    Sebastian Wallot \email{sebastian.wallot@ae.mpg.de} 
    }

\references{
  
  Webber Jr, C. L., and Zbilut, J. P. (2005). Recurrence
  quantification analysis of nonlinear dynamical systems.
  Tutorials in contemporary nonlinear methods for the
  behavioral sciences, 26-94.
  
  Marwan, N., and Kurths, J. Nonlinear analysis of bivariate
  data with cross recurrence plots. Physics Letters A 302.5
  (2002): 299-307.
  
  Coco, M. I., Monster, D., Leonardi, G., Dale, R., & Wallot, S. (2021). 
  Unidimensional and Multidimensional Methods for Recurrence 
  Quantification Analysis with crqa. R Journal, 13(1).
  
}

\keyword{ package }

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
           mindiagline, minvertline, tw, whiteline, recpt, side, method, metric, 
           datatype)

print(ans[1:10]) ## last argument of list is the recurrence plot


}
