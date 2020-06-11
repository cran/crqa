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
  
ts1 = as.vector(as.matrix(ts1));   ts2 = as.vector(as.matrix(ts2));
dpoint = seq(1, (length(ts1) - (windowsize)-1), windowstep)

drpd = vector() ## a vector to store all the points
i = 1
for (i in dpoint){
  
  ts1win = ts1[i:(i+windowsize - 1)];
  ts2win = ts2[i:(i+windowsize - 1)];
  
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
