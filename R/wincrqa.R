## wincrqa(): Function to run windowed recurrence analyses
##
## Arguments:
##      ts1                     : numeric vector or matrix ; first time series or full recurrence matrix
##      ts2                     : numeric vector ; second time series
##      windowstep              : numeric ; the interval by which the window is moved
##      windowsize              : numeric ; the size of the window
##      delay                   : numeric ; the lag for creating copies in phase-space reconstruction (i.e., lag intervals)
##      embed                   : numeric ; the number of embedding dimensions in phase-space reconstruction (i.e., number of copies)
##      radius                  : numeric ; the threshold distance by which points are considered recurrent;
##                                (default: .001)
##      rescale                 : numeric; choice for rescaling the distance matrix:
##                                    rescale = 0: no rescaling;
##                                    rescale = 1: rescale by mean distance of the entire matrix;
##                                    rescale = 2: rescale by maximum distance of the entire matrix;
##                                    rescale = 3: rescale by minimum distance of the entire matrix;
##                                    rescale = 4: rescale by Euclidean distance of the entire matrix;
##                                (default: 0)
##      normalize               : numeric ; choice for normalizing the time series:
##                                    normalize = 0: no normalizing;
##                                    normalize = 1: normalize by unit interval;
##                                    normalize = 2: normalize by z-score;
##                                (default: 0)
##      mindiagline             : numeric ; minimum number of sequential diagonal points to be considered a line
##                                (default: 2)
##      minvertline             : numeric ; minimum number of sequential vertical points to be considered a line
##                                (default: 2)
##      tw                      : numeric ; the Thieler window parameter, which removes the center diagonal line(s)
##                                (default: 0)
##      whiteline               : boolean ; whether or not to calculate empty vertical lines
##                                (default: FALSE)
##      recpt                   : boolean ; whether the input (ts1) is a recurrence matrix or not
##                                (default: FALSE)
##      side                    : character ; choice of side for the recurrence plot the measures are calculated (relative to center line):
##                                    side = "upper": upper half of the RP (from the center line and excluding the center line)
##                                    side = "lower": lower half of the RP (from the center line and excluding the center line)
##                                    side = "both": entire RP (including the center line)
##                                (default: "both")
##      method                  : character ; choice of recurrence analysis to perform:
##                                    method = "rqa": autorecurrence
##                                    method = "crqa": cross-recurrence
##                                    method = "mdcrqa": multidimensional cross-recurrence
##                                (default: "crqa")
##      metric                  : character ; choice of the type of distance metric used; argument will accept any method from rdist()
##                                (default: "euclidean")
##      datatype                : character ; the kind of data being analyzed:
##                                    datatype = "categorical": categorical data
##                                    datatype = "continuous": continuous data
##                                (default: "continuous")
##      trend                   : boolean ; whether or not to calculate the change in density of recurrence from the center line outward
##                                (default: FALSE)
## Value:
##     A dataframe of recurrence quantification metrics for each window
##
## Author(s): Moreno I. Coco 
## Contributed by: Alexandra Paxton
## Acknowledgements: based on the original MATLAB code by Rick Dale

.packageName <- 'crqa'

wincrqa <- function(ts1, 
                    ts2, 
                    windowstep, 
                    windowsize, 
                    delay, 
                    embed,
                    radius = 0.001, 
                    rescale = 0,
                    normalize = 0, 
                    mindiagline = 2, 
                    minvertline = 2,
                    tw = 0, 
                    whiteline = FALSE, 
                    recpt = FALSE, 
                    side = 'both', 
                    method = 'crqa',
                    metric = 'euclidean', 
                    datatype = 'continuous', 
                    trend = FALSE){
  
  ## stop immediately if the windowsize is smaller than delay*phase AND we're not supplying an RP
  if ((windowsize < embed*delay) & recpt == FALSE){ 
    stop("Phase-space (embed*delay) longer than windowsize")
  }  
  
  # print(data.frame(windowstep, windowsize, delay, embed, radius, rescale, 
  #                 normalize, mindiagline, minvertline, tw, whiteline, 
  #                 recpt, side, method, metric, datatype, trend))
  
  # initialize a flag that will be used to check whether the last window is included
  irregular = FALSE
  
  ## check the different contexts in which the analyses are running
  if (recpt == FALSE){
    
    # for RQA: check that the user as provided the same data to allow windowed recurrence
    if (method == "rqa"){
      if (exists("ts2")) {
        ts2 = ts2 
      } else {
        stop("Please provide the same vector as argument of t2")
      }
      maxd = length(ts1) ## the total number of points
    }
    
    # for CRQA
    if (method == "crqa"){
      ts1 = as.vector(as.matrix(ts1))
      ts2 = as.vector(as.matrix(ts2))
      maxd = length(ts1) ## the total number of points
      
    }
    
    # for MdCRQA
    if (method == "mdcrqa"){
      ts1 = as.matrix(ts1)
      ts2 = as.matrix(ts2)
      maxd = nrow(ts1) # the total number of points
    }
  } else {
    # for calculating from a supplied RP matrix
    logical_RP = as.logical(ts1)
    v1l = nrow(ts1)
   # print(logical_RP)
    ts1 = matrix(logical_RP, 
                 ncol = as.numeric(v1l))
    # print(dim(ts1))
    v1l = nrow(ts1)
    v2l = ncol(ts1)
    ind = which(ts1 > 0, arr.ind = TRUE)
    r = ind[,1]
    c = ind[,2]
    
    # total possible number of points on the RP
    maxd = v1l
  }
  
  
  ## need to make sure to include also "irregular" final segments
  points = seq(1, (maxd - (windowsize)-1), windowstep)
  
  fpoint = points[length(points)] ## beginning of last window
  
  ## as the dimension of the data may often be not that precise
  if (fpoint != (fpoint + windowsize)){
    irregular = TRUE ## flag used later on
    fpoint = fpoint + windowsize # update last point
    points = c(points,  fpoint) ## add final point
  }
  
  
  crqwin = vector() ## initialize the matrix to be filled with crq results
  tsp = 0 ## set a counter with all windows at which rec was computed    
  
  i = 1
  for (i in points){
    tsp = tsp +1
    
    if (irregular == TRUE & i == fpoint){ ## this is the last window to consider 
      ixs = i:maxd 
      # warning(paste("Your latest window was shorter: ", length(ixs)))
    } else {
      ixs = i:(i + windowsize - 1)
    }
    
    if (recpt == FALSE){
      if (method == "crqa" | method == "rqa"){
        ## the data consists of two vectors  
        ts1win = ts1[ixs];
        ts2win = ts2[ixs];
      }
      
      if (method == "mdcrqa"){
        ## the data consists of two matrices 
        ts1win = ts1[ixs, ];
        ts2win = ts2[ixs, ];
      }
    } else {
      # the data come from an already-constructed RP
      # print("RP inputted")
      ts1win = as.matrix(ts1[ixs, ixs])
      ts2win = NA
    }
    
    ## check whether the windowsize can accommodate the delay*embed OR the data come from a plot
    if ((length(ts1win) > delay*embed) || recpt == TRUE){
      
      # print(ts1win)
      
      ans = crqa(ts1win, ts2win, delay, embed, rescale, radius, 
                 normalize, mindiagline, minvertline,
                 tw, whiteline, recpt, side, method, metric,
                 datatype)
      
      RP = ans$RP
      
      # if we have an empty recurrence plot, do this
      if (length(RP) == 1) { RP = vector() }
      
      ## if trend needs to be calculated do it here
      ans = as.numeric( unlist(ans[1:10]) )
      
      #print(c(i, ans))
  
      if (trend == TRUE){
        
        if (length(RP) > 0){
          RP = as.matrix(RP) ## diags() has changed behaviour 
          NX = ncol(RP)
          T = vector("numeric", length = NX-1)
          
          ## below some test matrices
          # RP = matrix(1:16,nrow=4)
          # RP = cbind(c(2,7), c(1,9), c(10, 12))
          
          # do we need to span the full matrix or just half?
          ## create a vector of indeces for all diagonals of the matrix  
          # ixi_diag = -(ncol(RP)-1):(nrow(RP)-1) 
          
          # k = 2
          for (k in 1:(NX-1)){  ##ixi_diag
            # print(diags_extract(RP, k)) ## check the diagonal
            ixi_rec = which(diags_extract(RP, k) != F) ## find recurrence
            T[k] = length(ixi_rec)/ (NX-k)*100;
          }
          
          Ntau = NX - 1 - round(0.1*NX);
          
          ## last 10% of the RP will be skipped
          p = polyfit(2:(Ntau+1),T[1:Ntau], 1) # slope
          TREND = 1000 * p[1]
          ## Webber's definition includes factor 1000
          
        } else { TREND = NA}
      } else {
        TREND = NA
      }
      
      crqwin = rbind(crqwin, c(ans, tsp, TREND), deparse.level = 0)
      
    } else {
      warning( paste("Window", tsp, "was removed because windowsize < delay*embed", sep = " ") )
      
    }
  }
  
  ## name the measures
  colnames(crqwin) = c("RR", "DET", "NRLINE", "maxL", "L", 
                       "ENTR", "rENTR", "LAM", "TT", "catH", "win", "TREND")
  
  return(as.data.frame(crqwin))
  
}