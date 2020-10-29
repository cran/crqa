## written by Moreno I. Coco 2013 (moreno.cocoi@gmail.com)
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: step = interval considered on the serie;
##            windowstep = the step of the window. 
##            windowsize =  the size of the window
##           
##          
## other arguments to pass are the same as crqa

# windowsize = 200; windowstep = 50;
# type = 1; delay = 1; embed = 1;
# rescale = 1; radius = 0.0001;
# normalize = 0; minline = 2

# source("crqa.R")

# tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 1000)
# ts1 = tS[1,]; ts2 = tS[2,]

.packageName <- 'crqa'

wincrqa <- function(ts1, ts2, windowstep, windowsize, delay, embed,
                    radius = 0.001, rescale = 0,
                    normalize = 0, mindiagline = 2, minvertline = 2,
                    tw = 0, whiteline = FALSE, recpt = FALSE, side = 'both', 
                    method = 'crqa', metric = 'euclidean', 
                    datatype = 'continuous', trend = FALSE){
  
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
  
  crqwin = vector()
  tsp = 0 ## set a counter with the win at which rec was found    
  
  # i = 1621
  for (i in points){
    tsp = tsp +1
    
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
    
    ans = crqa(ts1win, ts2win, delay, embed, rescale,
               radius, normalize, mindiagline, minvertline,
               tw, whiteline, recpt, side, method, metric,
               datatype)
    
    RP = ans$RP
    if (length(RP) == 1) RP = vector() ## a trick for cases
    ## with empty recurrence plot
    ans = as.numeric( unlist(ans[1:10]) )
    ## if trend needs to be calculated do it here
    
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
        
        # k = 1
        for (k in 1:(NX-1)){  ##ixi_diag
#          print(diags(RP, k)) ## check the diagonal
          ixi_rec = which(diags(RP, k) != F) ## find recurrence
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
    
  }
  
  ## name the measures
  colnames(crqwin) = c("RR", "DET", "NRLINE", "maxL", "L", 
                       "ENTR", "rENTR", "LAM", "TT", "catH", "win", "TREND")
  
  return(as.data.frame(crqwin))
  
}

