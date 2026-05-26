\name{crqa-package}
\alias{crqa-package}

\docType{package}

\title{
  Unidimensional and Multidimensional Methods for Recurrence Quantification Analysis
}

\description{
  Auto, cross, and multidimensional recurrence quantification analysis
  (RQA/CRQA/MdCRQA). The package computes in-depth measures of recurrence
  plots for continuous and categorical time series, supports diagonal
  profile analysis, windowed analysis, and parameter optimisation.
  Since v2.1.0 the core kernel is implemented as a fused Rcpp inner
  loop with OpenMP parallelism (transparently falling back to serial
  on platforms without OpenMP support), and an approximative path
  (\code{method = "aRQA"}) scales to very long series (N >> 10^5)
  with negligible memory footprint.
  See Coco and Dale (2014) <doi:10.3389/fpsyg.2014.00510>,
  Wallot (2018) <doi:10.1080/00273171.2018.1512846>, and
  Coco et al. (2021) <doi:10.32614/RJ-2021-062> for methodological
  background.
}

\details{

\code{crqa}: Core recurrence function. Examines recurrent structures of a
  single (\code{method = "rqa"}), two (\code{method = "crqa"}), or
  multidimensional (\code{method = "mdcrqa"}) time series embedded in
  phase space. Since v2.1.0 also supports \code{method = "aRQA"} for
  approximative auto-RQA via phase-space histogram binning, scaling to
  N >> 10^5 in O(N) memory. A fused C++ inner loop (euclidean, maximum,
  manhattan metrics) avoids materialising the full distance matrix
  for the exact path as well.

\code{aRQA}: Approximative RQA via the phase-space histogram method of
  Schultz et al. (2015). Exported as a standalone function;
  also accessible via \code{crqa(\ldots, method = "aRQA")}.

\code{drpfromts}: Diagonal recurrence profile. Returns per-lag recurrence
  rates across a range of delays for auto, cross, or multidimensional RQA.

\code{lorenzattractor}: Simulate the Lorenz chaotic attractor.

\code{rosslerattractor}: Simulate the Roessler chaotic attractor.
  Particularly useful for benchmarking RQA implementations.

\code{mdDelay}: Estimate the time delay for embedding a multidimensional
  dataset using average mutual information.

\code{mdFnn}: Compute the percentage of false nearest neighbours for
  multidimensional time series as a function of embedding dimension.

\code{optimizeParam}: Iterative grid search over delay, embedding dimension,
  and radius to find parameters that yield a target recurrence rate.

\code{piecewiseRQA}: Block-wise RQA for structural decomposition of the RP.
  Splits the time series into consecutive blocks and aggregates the
  resulting measures. Useful when the research question calls for
  examining changes in recurrence structure across time, or when
  integrating RQA with sliding-window analyses at a coarser resolution
  than \code{wincrqa}. Supports parallel computation via \code{workers}.

\code{plot_rp}: Plot a recurrence plot returned by \code{crqa} using
  ggplot2.

\code{simts}: Simulate a coupled pair of categorical time series with
  controlled joint and conditional event probabilities.

\code{theiler_exclusion}: Helper to count the number of cells excluded
  by a Theiler window and one-sided mask, for use in the \code{rr_denom}
  calculation.

\code{wincrqa}: Windowed recurrence analysis. Computes a recurrence plot
  in overlapping windows and returns per-window recurrence measures.
  Supports parallel computation via \code{workers}.

\code{windowdrp}: Windowed diagonal recurrence profile. Like
  \code{wincrqa} but returns only the recurrence profile across delays
  averaged within each window.

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

print(ans[!names(ans) \%in\% "RP"])


}
