## plot the RP
## need to be improved, need to take out the for() loop for it
## RP  : the RP plot
## par : a list of parameters for the plotting
##     unit:   the unit to mark ticks on the axes
##     labelx: the label of the x-axis  
##     labely: the label of the y-axis  
##     cols  : the color for the recurrent recurrent points
##     pcex   : the size of the dots

.packageName <- 'crqa'

plotRP <- function(RP, par){
  
  if (exists("par") == FALSE){ # we use some defaults
    ## default values
    unit  = 2; labelx = "Time"; labely = "Time" 
    cols  = "black"; pcex = .3; pch = 1; las = 0;
    labax = seq(0, nrow(RP), unit); labay = seq(0, nrow(RP), unit);
  } else { # we load the values that we desire
    for (v in 1:length(par)) assign(names(par)[v], par[[v]])
  }
  
  xdim   = nrow(RP)
  ydim   = ncol(RP)
  
  RP = matrix(as.numeric(RP), nrow = xdim, ncol = ydim) # transform it for plotting
  
  ind = which(RP == 1, arr.ind = T)
  
  tstamp = seq(0, xdim, unit)
  
  par(mar = c(3.8, 3.8, 0.2,2), font.axis = 2, cex.axis = 1,
      font.lab = 2, cex.lab = 1.2)
  
  plot(tstamp, tstamp, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  matpoints(ind[,1], ind[,2],  cex = pcex, col = cols, pch = pch) 
  
  mtext(labelx, at = mean(tstamp), side = 1, line = 2.2, cex = 1.2, font = 2)
  mtext(labely, at = mean(tstamp), side = 2, line = 2.2, cex = 1.2, font = 2)
  
  
  #  if (is.numeric(labax)){ ## it means there is some default
  #    mtext(labax, at = seq(1, nrow(RP), nrow(RP)/10), side = 1, line = .5, cex = 1, font = 2)
  #    mtext(labay, at = seq(1, nrow(RP), nrow(RP)/10), side = 2, line = .5, cex = 1, font = 2)
  #  } else{
  mtext(labax, at = tstamp, side = 1, line = .5, cex = .8, font = 2, las = las)
  mtext(labay, at = tstamp, side = 2, line = .5, cex = .8, font = 2, las = las)
  
  # }
  
}
