\name{plot_rp}
\alias{plot_rp}

\title{
Plot a recurrence matrix 
}

\description{
Given a recurrence plot object (RP) returned from the function  \code{crqa()}, 
visualize it.
}

\usage{
plot_rp(rp_matrix,
        title = "",
        xlabel = "Time",
        ylabel = "Time",
        pcolour = "black",
        geom = c("tile", "point", "void"),
        flip_y = FALSE)
}

\arguments{
  \item{rp_matrix}{A recurrence plot sparse matrix from \code{crqa()}}
  \item{title}{Main title of the plot (character). Default: \code{""} (empty)}
  \item{xlabel}{The x-axis label (character). Default: \code{"Time"}}
  \item{ylabel}{The y-axis label (character). Default: \code{"Time"}}
  \item{pcolour}{The colour used for points in the plot (character). Default: \code{"black"}}
  \item{geom}{The ggplot2 geom used to draw the points (character).
  One of: \code{"tile"}, \code{"point"} or \code{"void"}.
  The default value, \code{"tile"}, uses (\code{geom_tile()} to plot 
  the points.
  The value \code{"point"} uses \code{geom_point()} to plot the points, and
  the value \code{"void"} produces a plot without any geom, so users can
  add their own for a more costumized plot (see examples below).}
  \item{flip_y}{Logical argument controlling whether to flip the directionality of
  the y axis (logical). Default: \code{FALSE}}
}

\details{
The argument \code{geom} is a character denoting what type of ggplot2 geom is
used to plot the points in the recurrence plot. Allowed values are 
\code{"tile", "point"} and \code{"void"}. The latter produces an empty plot, 
so the user is required to add their own geom to the output. See example below.
}

\value{
A ggplot2 plot.
}


\author{
Mønster, D. \email{danm@econ.au.dk}
}

\seealso{
This plotting solution as substituted what was previously \code{plotRP()} 
available till version crqa.2.0.5
}

\examples{
# Here is a very short example based on these two time series
ts_1 <-  c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1)
ts_2 <-  c(1, 1, 2, 2, 0, 0, 1, 2, 2, 1)

# Create the cross recurrence matrix
short_rec <-  crqa(ts_1, ts_2, 
                   delay = 1, embed = 1, radius = 0.001,
                   rescale = 1, method = "crqa")

# Extract the cross recurrence matrix
small_rp <-  short_rec$RP

# Make a cross recurrence plot with blue tiles
plot_rp(small_rp, pcolour  = "blue", geom = "tile")

# Make a cross recurrenbce plot with blue points
plot_rp(small_rp, pcolour  = "blue", geom = "point")

## Uncomment below if you want to have more example. Not runnable because
## of CRAN, globally unavailable function plus timing constraint to check vignette.

# Neither of the two plots above may be suitable
# Using geom = "void" we can add a custumized geom, here blue tiles
# that are smaller than the default width and height
# plot_rp(small_rp, geom = "void") +
#  geom_tile(fill = "blue", height = 0.75, width = 0.75)

# Another custom geom uses tiles with a border colour
# plot_rp(small_rp, geom = "void") +
#  geom_tile(fill = "blue", colour = "pink", linewidth = 1)

# Use the eyemovement dataset to create a cross-recurrence plot
# large_crqa <- crqa(eyemovement$narrator, eyemovement$listener,
#                delay = 1, embed = 1, radius = 0.001,
#                method = "crqa", metric = "euclidean", datatype = "continuous")

# Extract the recurrence matrix from the result
# large_rp <- large_crqa$RP

# Create a recurrence plot using defaults
# plot_rp(large_rp)

# The same recurrence plot with flipped y axis and using geom_point()
# plot_rp(large_rp, flip_y = TRUE, geom = "point", title = "Flipped y axis") 

# Add axes labels, a title and extra plot elements
# Note that we can add to the output, here added the line of synchonization
# using ggplot2's geom_abline() function.
# plot_rp(large_rp, 
#        xlabel = "Narrator",
#        ylabel = "Listener",
#        title = "Coloured tiles with line of synchronization",
#        pcolour = "blue") +
#  geom_abline(intercept = 0, slope = 1)

}
