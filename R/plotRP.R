## plot the RP

## RP  : the RP plot
## par : a list of parameters for the plotting
##     unit:   the unit to mark ticks on the axes
##     labelx: the label of the x-axis  
##     labely: the label of the y-axis  
##     cols  : the color for the recurrent recurrent points
##     pcex   : the size of the dots

.packageName <- 'crqa'

plotRP <- function(RP, par){
  
  unit = cols = labelx = labely = pcex= NA
  
  for (v in 1:length(par)) assign(names(par)[v], par[[v]])
  
  xdim   = nrow(RP)
  ydim   = ncol(RP)
  
  RP = matrix(as.numeric(RP), nrow = xdim, ncol = ydim) # transform it for plotting
  
  tstamp = seq(0, xdim, unit)
  
  par(mar = c(3.8, 3.8, 0.2,2), font.axis = 2, cex.axis = 1,
      font.lab = 2, cex.lab = 1.2)
  
  plot(tstamp, tstamp, type = "n", xlab = "", ylab = "")
  
  l = 1
  for (l in 1:ydim){
    ind = which(RP[,l] == 1)
    points(rep(l,length(ind)), ind, cex = pcex, col = cols, pch = 20)
  }
  
  mtext(labelx, at = mean(tstamp), side = 1, line = 2.2, cex = 1.2, font = 2)
  mtext(labely, at = mean(tstamp), side = 2, line = 2.2, cex = 1.2, font = 2)
  
}
