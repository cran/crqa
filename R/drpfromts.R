# GNU License: written by Moreno I. Coco (moreno.cocoi@gmail.com)
# based on code of Rick Dale (rdale@ucmerced.edu)

# compute diagonal recurrence profile for two
# categorical time series, t1 t2 (x / y)
# the windowsize is the total nr. of lags over which the series are evaluated

# modified as of crqa.1.0.7, now we use the crqa function and directly
# extract from it the diagonal means. This change was implemented to obtain 
# diagonal profiles from delayed and embedded continuous time-series

# the drpfromts assumes as default categorical variables, 
# no-delay, no-embed, and radius < 0.001

.packageName <- 'crqa'

drpfromts <- function(ts1, ts2, windowsize,  
                      radius = 0.001, delay = 1, embed = 1, rescale = 0,
                      normalize = 0, mindiagline = 2, minvertline = 2,
                      tw = 0, whiteline = F, recpt = F, side = 'both', 
                      method = 'crqa', metric = 'euclidean', 
                      datatype = 'categorical'){
  
  
  if(method != "mdcrqa"){ ## in case the method is multi-dimensional do not perform checks on vectors
    
    if (datatype == "categorical"){ 
      ## convert the data into numeric
      ts1     = as.character( as.matrix( ts1 ) )
      ts2     = as.character( as.matrix( ts2 ) )
      
      ## apply convenient function to transform the levels of the categorical variables into numerics
      tsnorm = numerify(ts1, ts2)
      ts1 = tsnorm$nwts1
      ts2 = tsnorm$nwts2
      
    }
    
    if (datatype == "continuous"){
      ## just double check that the data is numeric
      ts1 = as.numeric( as.matrix( ts1 ) )
      ts2 = as.numeric( as.matrix( ts2 ) )
      
    }
  }
  
  res = crqa(ts1, ts2, delay, embed, rescale, radius, normalize, 
             mindiagline, minvertline, tw, whiteline, recpt, side, 
             method, metric, datatype)
  
  RP = res$RP
  
  if (is.logical(RP) !=  T){ # we have some point that recur
    len      = nrow(RP)-1 # the nrow of the RP
    lags     = -len:len   # the diagonal of the RP
    
    RP      = matrix(as.numeric(RP), nrow = nrow(RP), ncol = ncol(RP))
    RP_form = unlist(lapply(split(RP, row(RP) - col(RP)), mean, na.rm = TRUE))
    wn      = which(lags >= -windowsize &  lags <= windowsize)
    
    drpd    = RP_form[wn]
  } else {
    drpd    = 0 # or technically it should be an NA
  }
  
  ## extract  max recurrence and the lag at which it occurred
  
  maxrec = max(drpd);
  maxlag = which(drpd == maxrec);
  
  return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )
  
}
