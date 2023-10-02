## plotRP(): Function to plot the recurrence plot
##
## Arguments:
##     RP          : the RP ngCMatrix output from crqa() function
##     par         : a list of parameters for the plotting, including:
##       unit      : numeric ; gap between sample labeling on axes.
##                   Note: only relevant if `show_ticks = TRUE`.
##                   (default: 10)
##       labelmain : character ; main title text of the plot 
##                   (default: "Recurrence Plot")
##       labelx    : character ; the text label of the x-axis 
##                   (default: "Time")
##       labely    : character ; the text label of the y-axis 
##                   (default: "Time")
##       cols      : character ; the color for the recurrent points;
##                   may include any colors from the base R plot repertoire
##                   (default: "black)
##       pcex      : numeric ; the size of the recurrent points
##                   (default: .3)
##       pch       : numeric ; the style of the recurrent points
##                   (default: 1)
##       show_ticks: boolean ; whether to show x- and y-ticks or not
##                   (default: FALSE)
##
## Value:
##     A square plot visualizing the recurrence matrix.
##
## Author(s): Moreno I. Coco & Alexandra Paxton

.packageName <- 'crqa'

plotRP <- function(RP, par){
  
  ## default values for plotting
  # initialize variable to go around global visible binding
  # specify the space of possible parameters
  
  labelx = "Time"; labely = "Time"; labelmain = "Recurrence Plot"; cols = "black" 
  pcex = .3; pch = .3; unit = 10; show_ticks = FALSE 
  
  default_par = data.frame(labelx, labely, labelmain, cols,
                  pcex, pch, unit, show_ticks)
  
  # if no user-defined parameters exist, use defaults
  if (exists("par") == FALSE){
    for (default_v in 1:length(default_par)){ 
      assign(names(default_par)[default_v], default_par[[default_v]])
    }
  } else {
    
    # define user-defined variables
    for (v in 1:length(par)){ assign(names(par)[v], par[[v]]) }
    
    # find which labels have not been user-defined and use defaults
    undefined_par = names(default_par)[!names(default_par) %in% names(par)]
    use_defaults = default_par[undefined_par]
    for (default_v in 1:length(use_defaults)){
      assign(names(use_defaults)[default_v], use_defaults[[default_v]])
      }
  }
  
  # find the size of the RP
  xdim   = nrow(RP)
  ydim   = ncol(RP)
  
  # transform the RP into a square matrix
  RP = matrix(as.numeric(RP), nrow = xdim, ncol = ydim)
  
  # figure out where the recurrent points are
  ind = which(RP == 1, arr.ind = T)
  
  # create all time series of all possible samples on axis
  tstamp = seq(0, xdim, 1)
  
  # create the shell of the plot
  par(mar = c(3.8, 3.8, 2,2), font.axis = 2, cex.axis = 1,
      font.lab = 2, cex.lab = 1.2)
  plot(tstamp, tstamp, type = "n", 
       xlab = "", ylab = "", main=labelmain,
       xaxt = "n", yaxt = "n",
       font = 2)
  
  # add recurrent points to the plot
  matpoints(ind[,1], ind[,2], cex = as.numeric(pcex), col = cols, pch = as.numeric(pch))
  
  # add x- and y-axis labels to the plot (so they're not too far away)
  mtext(labelx, at = mean(tstamp), side = 1, line = 1.1, cex = 1.2, font = 2)
  mtext(labely, at = mean(tstamp), side = 2, line = 1.1, cex = 1.2, font = 2)
  
  # if the user would like x- and y-tickmarks, show and label them
  if (show_ticks == TRUE){
    labax = seq(0, nrow(RP), unit)
    labay = seq(0, nrow(RP), unit)
    las = 0
    mtext(labax, at = labax, side = 1, line = .1, cex = .8, font = 2, las = las)
    mtext(labay, at = labay, side = 2, line = .1, cex = .8, font = 2, las = las)
  }
  
  # create variable to save plot
  output_plot = recordPlot()
  invisible(dev.off())
  
  # return the plot
  return(output_plot)
}
