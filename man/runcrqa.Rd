\name{runcrqa}
\alias{runcrqa}

\title{Wrapper to compute different types of
  cross-recurrence quantification analyses
  between two time-series.}

\description{
  
  Wrapper to extract CRQ information on categorical and
  continuous time series information. The function
  provides two types of analysis: the recurrence diagonal
  profile (type = 1), or a detailed analysis of the
  recurrence plot (type = 2). For both methods, the
  function can perform  a profile analysis
  (method = 'profile'), which looks at how recurrence
  change for the different lags, or a window analysis
  (method = 'window'), where recurrence is tracked across
  the time-course by sliding overlapping windows.

}
\usage{
runcrqa(ts1, ts2, par)
}

\arguments{

  \item{ts1}{First time-series}
  \item{ts2}{Second time-series}
  \item{par}{A list of argument parameters depending on
    whether the wrapper is used to obtain only profiles,
    i.e., type 1, or to extract more detailed measures
    from the cross-recurrence plot, type 2. See details
    below for a detailed explanation of the
    arguments for the two different methods.}

}

\details{
  
  Independently of the type, the argument 'method'
  can take two values either 'profile' or 'window':
  
  \code{method = 'profile'}: compute the recurrence profile
  over the all time series for different lags
  \code{method = 'window'}: compute recurrence over time by
  sliding a window
  
  For \code{type 1}:
  
  with \code{method = 'profile'}

  \code{ws} = the width (+/-) timestamps to use to lag
  the series

  with \code{method = 'window'}  

  \code{step} = the sampling jumps over which windows
  are rolled

  \code{windowsize} = the size of the window of analysis.

  \code{lagwidth} = the number of lags to be analyzed.

  For \code{type 2}:

  \code{delay} = the delay introduced to the time series

  \code{embed} = the embedding dimensions for phase-space
  reconstruction

  \code{rescale} = Normalization of the distance matrix
  if rescale = 1 (mean distance of entire matrix);
  if rescale = 2 (maximum distance of entire matrix)
  
  \code{radius} = Maximum distance to accept two-points
  as recurrent. If the series are categorical, it must be
  set to a very small value

  \code{normalize} = Rescale factor for time-series;
  if normalize = 0 (do nothing);
  if normalize = 1 (Unit interval);
  if normalize =  2 (z-score)
  
  \code{mindiagline} = A minimum diagonal length of
    recurrent points. Usually set to 2, as it takes
    a minimum of two points to define any line.
  
  \code{minvertline} = A minimum vertical length of
  recurrent points.
  
  \code{whiteline} = A logical flag to calculate (TRUE)
  or not (FALSE) empty vertical lines.
  
  \code{recpt} = A logical flag indicating whether
    measures of cross-recurrence are calculated directly
    from a recurrent plot (TRUE) or not (FALSE).
    
  }

\value{

  The values returned depends on the type (1,2) of
  computation requested.

  For \code{type = 1}:

  \code{profile} = A vector of recurrence
  (ranging from 0,1) with length equal
  to the number of delays explored

  \code{maxrec} = Maximal recurrence observed between
  the two-series
  \code{maxlag} = Delay at which maximal recurrence
  is observed

  For \code{type = 2}:

  \code{rec} = The percentage of recurrent points
  falling within the specified radius. 

  \code{det} = Proportion of recurrent points forming
  diagonal line structures.
  
  \code{nrline} = The total number of lines in the
  recurrent plot

  \code{maxline} = The length of the longest diagonal
  line segment in the plot, excluding the main diagonal

  \code{meanline} = The average length of line structures

  \code{entropy} = Shannon information entropy of all
  diagonal line lengths
  
  \code{relEntropy} = Entropy measure normalized by the
  number of lines observed in the plot.
  Handy to compare across contexts and conditions

  \code{lam} = Proportion of recurrent points forming
  vertical line structures
  
  \code{tt} = The average length of vertical line structures

  
  }
  
\references{
  
  Webber Jr, C. L., and Zbilut, J. P. (2005). Recurrence
  quantification analysis of nonlinear dynamical systems.
  Tutorials in contemporary nonlinear methods for the
  behavioral sciences, 26-94.

  Marwan, N., and Kurths, J. Nonlinear analysis of
  bivariate data with cross recurrence plots.
  Physics Letters A 302.5 (2002): 299-307.

}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}


\seealso{ \code{\link{drpdfromts}}, \code{\link{windowdrp}},
\code{\link{crqa}}, \code{\link{wincrqa}}}

\examples{

data(crqa)

####################################################
###Cross-recurrence diagonal profile 

par = list(type = 1, ws = 100, method = "profile",
    datatype = "categorical", thrshd = 8, radius = 2,
     pad = FALSE)

ans = runcrqa(RDts1, RDts2, par)

profile = ans$profile; maxrec = ans$maxrec; maxlag = ans$maxlag

####################################################
###Windowed cross-recurrence profile 

par = list(type = 1, windowstep =  20, windowsize =  100, 
  method = "window", datatype = "categorical", thrshd = 8,
     pad = FALSE)

ans = runcrqa(RDts1, RDts2, par)

print(ans)

####################################################
### Cross-recurrence measures 

par = list(type = 2, delay = 1, embed =  1, rescale =  1,
  radius = 0.00001, normalize = 0, mindiagline = 2,
   minvertline = 2, whiteline = FALSE, recpt = FALSE)

res = runcrqa(RDts1, RDts2, par)

res[1:9]

}
