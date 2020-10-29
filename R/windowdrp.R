## GNU License, written by Moreno I. Coco
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: windowstep = interval considered on the serie;
##            window_size = the size of the window wherin crqa is runned. 
##            lag_profile_width = lags within the window

## windowstep = 10;
## windowsize = 100;
## lagwidth = 20;
## datatype = "categorical"

.packageName <- 'crqa'

windowdrp <- function(ts1, ts2, windowstep, windowsize, lagwidth, 
                      radius = 0.001, delay = 1, embed = 1, rescale = 0,
                      normalize = 0, mindiagline = 2, minvertline = 2,
                      tw = 0, whiteline = FALSE, recpt = FALSE, side = 'both', 
                      method = 'crqa', metric = 'euclidean', 
                      datatype = 'continuous'){
  
  irregular = FALSE; # initialize a flag that will be used to check whether the last window is included
  ## check the different contexts in which the analyses are running
  
  if (method == "rqa"){
    ## check that the user as provided the same data to allow windowed recurrence
    if (exists("ts2")) ts2 = ts2 else stop("Please provide the same vector as argument of t2")
  }
  
  if (method == "crqa"){
    ts1 = as.vector(as.matrix(ts1));   ts2 = as.vector(as.matrix(ts2))
    maxd = length(ts1) ## the total number of points
    points = seq(1, (maxd - (windowsize)-1), windowstep)
  }
  
  if (method == "mdcrqa"){
    ts1 = as.matrix(ts1);   ts2 = as.matrix(ts2)
    maxd = nrow(ts1) ## the total number of points
    points = seq(1, (maxd - (windowsize)-1), windowstep)
  }
  
  ## need to make sure to include also "irregular" final segments
  fpoint = points[length(points)] ## beginning of last window
  ## as the dimension of the data often may not be that precise
  if (fpoint != (fpoint + windowsize)){
    irregular = TRUE ## flag used later on
    fpoint = fpoint + windowsize # update last point
    points = c(points,  fpoint) ## add final point
  }
  
  
  # i = 1
  drpd = vector() ## a vector to store all the points
  
  for (i in points){
    if (irregular == TRUE & i == fpoint){ ## this is the last window to consider 
      ixs = i:maxd 
      # warning(paste("Your latest window was shorter: ", length(ixs)))
    } else {
      ixs = i:(i + windowsize - 1)
    }
    
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
    
    drpd = c(drpd,
             mean(drpfromts(ts1win, ts2win,lagwidth,
                            radius, delay, embed, rescale,
                            normalize, mindiagline, minvertline,
                            tw, whiteline, recpt, side, 
                            method, metric, 
                            datatype)$profile))
    
  }
  
  maxrec = max(drpd);
  maxlag = which(drpd == maxrec);
  
  return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )
  
}
