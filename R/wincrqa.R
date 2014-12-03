## GNU License, written by Moreno I. Coco
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: step = interval considered on the serie;
##            window_size = the size of the window wherin crqa is runned. 
##            lag_profile_width = lags within the window
##            opt = list( type = 1 - do, or 0 - don't do crqa on windows;
##            if type 1; opt should also specify parameters to crqa 
##            (delay, embed, rescale, radius, normalize, minline)
##            
# step = 10;
# windowsize = 100;
# lagwidth = 20;
# type = 1; delay = 1; embed = 1; rescale = 1; radius = 0.00001; normalize = 0; minline = 2

# source("crqa.R")

.packageName <- 'crqa'

wincrqa <- function(x, y, step, windowsize, lagwidth, delay, embed, 
rescale, radius, normalize, mindiagline, minvertline, tw, whiteline, recpt){

    recpt = FALSE
    ## we do not expect as input a windowed recurrent plot
    x = as.vector(as.matrix(x));   y = as.vector(as.matrix(y))
    points = seq(1, (length(x) - (windowsize)-1), step)
    
    crawin = vector()
    tsp = 0 ## set a counter with the win at which rec was found
    
    
    for (i in points){
        tsp = tsp +1
        
        xwin = x[i:(i+windowsize)];
        ywin = y[i:(i+windowsize)];
        
        ans = crqa(xwin, ywin, delay, embed, rescale,
            radius, normalize, mindiagline, minvertline, tw,
            whiteline, recpt)
        
    
        ans = as.numeric( unlist(ans[1:9]) )
        ## for window recurrent do not save the recurrence plot
        
        crawin = rbind(crawin, c(ans, tsp), deparse.level = 0)
        
    }
    
    
    return(crawin)
    
}

