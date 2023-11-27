\name{plotRP}
\alias{plotRP}

\title{Plot a recurrence matrix}

\description{
  A convenience function to plot the RP matrix returned by
  `crqa()`
}

\usage{ plotRP(RP, labelx, labely, labelmain, 
                   cols, pcex, pch, unit, 
                   show_ticks) }

\arguments{
  \item{RP}{A recurrence plot sparse matrix from `crqa()`}
  \item{labelx}{The text label of the x-axis (character). Default: "Time"}
  \item{labely}{The text label of the y-axis (character). Default: "Time"}
  \item{labelmain}{Main title text of the plot (character). Default: "Recurrence Plot"}
  \item{cols}{The color for the recurrent points (character). May include any colors from the base R plot repertoire. Default: "black"}
  \item{pcex}{The size of the recurrent points (numeric). Default: .3}
  \item{pch}{The style of the recurrent points (numeric). Default: 1}
  \item{unit}{Gap between sample labeling on axes (numeric). Note: only relevant if `show_ticks = TRUE`. Default: 10}
  \item{show_ticks}{Whether to show x- and y-ticks or not (boolean). Default: FALSE}
  }


\value{A square plot visualising the recurrence matrix, saved as a `output_plot` object
}

\author{Moreno I. Coco \email{moreno.cocoi@gmail.com}
Alexandra Paxton \email{alexandra.paxton@uconn.edu}
}

\examples{
## run a simple crqa
ts1 = c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1)
ts2 = c(1, 1, 2, 2, 0, 0, 1, 2, 2, 1)
delay = 1; embed = 1; rescale = 1; radius = 0.001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0;

res = crqa(ts2, ts1, delay, embed, rescale, radius,  normalize, mindiagline,
minvertline, tw)

RP = res$RP
plotRP(RP, labelx = "Time", labely = "Time", show_ticks = TRUE, unit = 2,
       cols = "blue", pcex = 1)
}

\keyword{misc}
