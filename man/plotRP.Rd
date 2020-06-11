\name{plotRP}
\alias{plotRP}

\title{Plot a recurrence matrix}

\description{
  A convenience function to plot the RP matrix returned by 
  the `crqa()'  
}

\usage{ plotRP(RP, par) }

\arguments{
  \item{RP}{The Recurrence Plot sparse matrix data}
  \item{par}{a list of parameters for the plotting:
            unit: to unit to define ticks on the axes relative to RP dimension 
            labelx: the label of the x-axis  
            labely: the label of the y-axis  
            cols: the color for the recurrent recurrent points
            pcex: the size of the dots
            pch: the point type
            labax = labels for the x-axis
            labay = labels for the y-axis
            las = orientation of labels}
            }

\value{
  A square plot visualising the recurrence matrix.
}

\author{Moreno I. Coco (moreno.cocoi@gmail.com)}


\examples{
## run a simple crqa 
ts1 = c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1)
ts2 = c(1, 1, 2, 2, 0, 0, 1, 2, 0, 0)

res = crqa(ts2, ts1, delay = 1, embed = 1, rescale = 0, 
           radius = 0.01,  normalize = 0, mindiagline = 2, 
           minvertline = 2, tw = 0, method = "crqa", side = "both", 
           datatype = "continuous")

RP = res$RP

## define plotting arguments
par = list(unit = 1, labelx = "Time", labely = "Time", 
           cols = "blue", pcex = 1, pch = 19, 
           labax = seq(0, nrow(RP), 1), 
           labay = seq(0, nrow(RP), 1), 
           las = 1)

plotRP(RP, par)

}

\keyword{misc}
