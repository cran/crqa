# GNU License: written by Moreno I. Coco (moreno.cocoi@gmail.com)
# based on code of Rick Dale (rdale@ucmerced.edu)

# compute diagonal recurrence profile for two
# categorical time series, t1 t2 (x / y)
# the ws is the total nr. of lags over which the series are evaluated

# modified as of crqa.1.0.7, now we use the crqa function and directly
# extract from it the diagonal means. This change was implemented to obtain 
# diagonal profiles from delayed and embedded continuous time-series

# the drpdfromts assumes as default categorical variables, 
# no-delay, no-embed, and radius < 0.001

.packageName <- 'crqa'

drpdfromts <- function(ts1, ts2, ws, datatype, 
                       radius = 0.001, delay = 1, embed = 1, rescale = 0,
                       normalize = 0, mindiagline = 2, minvertline = 2,
                       tw = 0){
  
  if (datatype == "categorical"){ 
    ## convert the data into numeric
    ts1     = as.character( as.matrix( ts1 ) )
    ts2     = as.character( as.matrix( ts2 ) )
    
    checked = checkts(ts1, ts2, datatype, thrshd = length(ts1), pad = FALSE)
    ts1     = checked[[1]][, 1]
    ts2     = checked[[1]][, 2]
    
  }
  
  if (datatype == "continuous"){
    ## just double check that the data is numeric
    ts1 = as.numeric( as.matrix( ts1 ) )
    ts2 = as.numeric( as.matrix( ts2 ) )
    
  }
  
  ## run full blown crqa
  res = crqa(ts2, ts1, delay, embed, rescale, radius, 
             normalize, mindiagline, minvertline, tw)
  
  ## extract the RP
  RP = res$RP
  if (class(RP) != "logical"){ # we have some point that recur
    len      = nrow(RP)-1 # the nrow of the RP
    lags     = -len:len   # the diagonal of the RP
    
    RP      = matrix(as.numeric(RP), nrow = length(ts1), ncol = length(ts2))
    RP_form = unlist(lapply(split(RP, row(RP) - col(RP)), mean, na.rm = TRUE))
    wn      = which(lags >= -ws &  lags <= ws)
    
    drpd    = RP_form[wn]
  } else {
    drpd    = 0 # or technically it should be an NA
  }
  
  ## extract  max recurrence and the lag at which it occurred
  
  maxrec = max(drpd);
  maxlag = which(drpd == maxrec);
  
  return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )
  
}

# rm(list = ls())

#data(crqa)

#ts1 = as.numeric(as.matrix(RDts1))
#ts2 = as.numeric(as.matrix(RDts2))

#res = drpdfromts(ts1, ts2, ws = 100, datatype = "continuous", radius = 0.001, 
#                 delay = 1, embed = 1, rescale = 0, normalize = 0, 
#                 mindiagline = 2, minvertline = 2, tw = 0)

#profile = res$profile

#plot(seq(1,length(profile),1), profile, type = "l", lwd = 5,
#     xaxt = "n", xlab = "Lag", ylab = "Recurrence")


