# A ggplot2-based function to produce a recurrence plot from a
# matrix. Much faster than using graphics::image() and nicer than
# crqa::plotRP().
# Author Dan Moester
.packageName <- 'crqa'

plot_rp <- function(rp_matrix,
                    title = "",
                    xlabel = "Time",
                    ylabel = "Time",
                    pcolour = "black",
                    geom = c("tile", "point", "void"),
                    flip_y = FALSE) {
  
  ## make globally visible some objects that will later on be filled in
  Var1 = Var2 = value = NA
  
  geom <- match.arg(geom)
  # In most cases rp_matrix is a sparse matrix from the Matrix package,
  # so we convert it to a regular matrix first.
  rp_matrix <- as.matrix(rp_matrix)
  # Get the dimensions for check and to use in the plot
  rp_dims <- dim(rp_matrix)
  if (rp_dims[1] != rp_dims[2]) {
    stop("Matrix is not square. Number of rows should be equal to number of columns.")
  }
  
  # Create a data frame, that can be passed to ggplot2
  #
  # Only keep recurrent points for better performance, since the non-recurrent
  # points are not plotted anyway.
  #
  # The row and column indices are returned as factors by the function
  # as.data.frame.table, so they have to be re-coded as integers.
  #
  # Finally the recurrent points are changed from numeric (1) to a factor
  # to get a discreet variable.
  #
  # NOTE: To be compatible with R version 3, no native pipes are used, but
  # repeated function calls are used instead.
  rp_df <- as.data.frame.table(rp_matrix, responseName = "value")
  rp_df <- filter(rp_df, value == 1)
  rp_df <- mutate(rp_df, 
                  Var1 = as.integer(Var1),
                  Var2 = as.integer(Var2),
                  value = factor(1L))
  # Create a rectangular plot and set the limits. Geom will be added later
  rp <- ggplot(rp_df, aes(x = Var1, y = Var2, fill = value)) +
    coord_fixed(ratio = 1) +
    xlab(xlabel) +
    ylab(ylabel) +
    labs(title = title) +
    expand_limits(x = 1, y = 1) +
    expand_limits(x = rp_dims[1], y = rp_dims[2]) +
    theme_classic() +
    theme(legend.position = "none")
  # Flip the y axis if requested
  if (flip_y) {
    rp <- rp + scale_y_reverse()
  }
  # Now add the geom, either tile or point
  if (geom == "tile") {
    rp <- rp + 
      scale_fill_manual(values = c("0" = "white", "1" = pcolour)) +
      geom_tile(width = 0.9, height = 0.9) 
  } else if (geom == "point") {
    rp <- rp +
      scale_color_manual(values = c("0" = "white", "1" = pcolour)) +
      geom_point(aes(colour = value),
                 size = 0.5, stroke = 0, shape = ".")
  } else if (geom == "void") {
    
  }
  rp
}
