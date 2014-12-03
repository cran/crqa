# GNU > 2 License:
## written by Moreno I. Coco (moreno.cocoi@gmail.com)

# Iterative procedure exploring a combination of parameter values
# to obtain maximal recurrence between two time-series
# The procedure should be applied on several trials to find the
# "optimal" setting for the parameters...
# it should be applied on time-series with continuous data

# TODO: probably, for windows the selection of size/embedding
# has to be done at the same time. 


## source("./ami.R")       ## average mutual information
## between 2 series
## source("./crqa.R")

## lgM = 20 ## maximum lag to check for average mutual information
## steps = seq(1, 10, 1) ##how many points we should look
## ahead to define a local minimum
  
## cut.del = seq(1, 40,1)  ## cut off of delay is set
## radiusspan = 2 ## increasing radius step means smaller steps
## radiussample = 10 ## number of radius points within the steps
                     ## to be sampled

.packageName <- 'crqa'

optimizeParam <- function(ts1, ts2, par){
    rescale = normalize = mindiagline = minvertline =
        lgM = steps = cut.del = recpt =  whiteline = tw =
            radiusspan = radiussample = NULL
    ## initialize attached variables
    
    for (v in 1:length(par)) assign(names(par)[v], par[[v]]) ## assign parameters
    
    if (radiusspan <= 1){
        stop("Radius span too small, please choose a larger value")
    }
                                        #  require(tseriesChaos)
    
    ## check that the input is vector and numeric
    
    if (is.vector(ts1) != TRUE){ ts1 = as.numeric(as.matrix(ts1) ) }
    if (is.vector(ts2) != TRUE){ ts2 = as.numeric(as.matrix(ts2) ) }

    
    mi1 = as.numeric( mutual(ts1, lag.max = lgM, plot = FALSE))
    mi2 = as.numeric( mutual(ts2, lag.max = lgM, plot = FALSE))
    
    ## iterate over a range of possible steps to find out the optimal one
    
    lag1 = lag2 = vector()
    for (s in steps){
                                        # print(s)
        
        for (l in 1:(lgM - s)){
            if (mi1[l] < mi1[l + s]){
                lg1 = which(mi1 == mi1[l])
                lag1 = c(lag1, lg1)
        break  }
        }
        
        for (l in 1:(lgM - s)){
            if (mi2[l] < mi2[l + s]){
                lg2 = which(mi2 == mi2[l])
                lag2 = c(lag2, lg2)
                break  }
        }
    }
    
    if (length(lag1) == 0 | length(lag2) == 0){
        stop("Please try varying maximum lag (and/or step size):
      minimal mutual information was not reached")
    }
    
    ## take unique lags, i.e., avoid repeat
    lag1 = unique(lag1); lag2 = unique(lag2)
    
    
    ## for each time-series take the lag where MI is minimal
    
    ix.lg1 = which( mi1[lag1] == min(mi1[lag1]) )
    if (length(ix.lg1) > 1){
        ix.lg1 = sample(ix.lg1, 1) ## if there are more than one minimum
                                 ## sample one 
        lg1 = lag1[ix.lg1]}
    
    ix.lg2 = which( mi2[lag2] == min(mi2[lag2]) )
    if (length(ix.lg2) > 1){
        ix.lg2 = sample(ix.lg2, 1)
        lg2 = lag2[ix.lg2]}
    
    
#  mi = ami(ts1, ts2, lag = 1)[[1]] ## not sure about this function
#                                   ## perhaps look at it later

#  Choose a delay that will accommodate both ts.
#  For example if one has a considerably longer delay
#  indicated than another, you may need to select the longer
#  delay of the two because that ensures that new information
#  is gained for both by using the higher delay.
#  If delays are close to each other, you may want to use a
#  delay somewhere in between the two.
          
    delay.in = vector()
    
    for (delay in cut.del){
        ## print(delay)
        
        if (lg1 > lg2 & abs(lg1 - lg2) > delay){
          ## print( c(delay, lg1) )  
            delay.in = c(delay.in, lg1)
        }
        
        else if (lg2 > lg1 & abs(lg2 - lg1) > delay){
            ## print( c(delay, lg2) )            
            delay.in = c(delay.in, lg2)
        }
        
        else {
            delay.in = c(delay.in,
                round( mean(c(lg1, lg2))) )
        }
        
        ## print(delay.in)
        
    }
    
    ## take the delay that more frequently is found
    ## to be optimal
    
    del.in = as.data.frame( table(delay.in) )
    del = as.numeric( as.vector(
        del.in[which(del.in[,2] == max(del.in[,2])),1]))
    
    if (length(del) > 1){ del = del[length(del)] } ## take longer delay 
    
    
    ## Determine embedding dimension: %FNN function bottoms out
    ## (i.e., where it buys you nothing to add more dimensions).
    ## If the embedding dimension indicated for the two files is
    ## different, you’ll want to go with the higher embedding
    ## dimension of the two to make sure that you have
    ## sufficiently unfolded both time series.
    ## It’s generally safe to overestimate embedding dimension.
    ## Though this can amplify the influence of noise for
    ## high embedding dimensions. 
    
    ## TODO: bring out this parameters to optimize them also    
# m	 maximum embedding dimension = 20
# d	 delay parameter = get it from previous optimization
# t	 Theiler window = 0 ## in cross-recurrence LOC is important
# rt	 escape factor = leave it default
# eps	 neighborhood diameter = leave default
            
    embdts1 = false.nearest(ts1, m = 20, d = del, t = 0,
        rt = 10, eps = sd(ts1)/10)
    embmints1 = as.numeric( which(embdts1[1,] == min(embdts1[1,],
                                     na.rm = TRUE) ) )
    
    embdts2 = false.nearest(ts2, m = 20, d = del, t = 0,
        rt=10, eps=sd(ts2)/10)
    embmints2 = as.numeric( which(embdts2[1,] == min(embdts2[1,], na.rm = TRUE) ) )
    
    if ( length(embmints1) > 1){ embmints1 = embmints1[1] }
    if ( length(embmints2) > 1){ embmints2 = embmints2[1] }
    
    embdim = max( c(embmints1, embmints2) )        
    
    ## take the optimal parameters from above and then check
    ## set the radius
    
    ## first we need to see in the rescaled matrix what is
    ## the range of values, so that the radius interval
    ## can be consciously estimated.

    dm = rdist(ts1, ts2)
    if (rescale > 0){
            switch(rescale,
                   {1  ## Create a distance matrix that is re-scaled
                    ## to the mean distance
                    
                    rescaledist = mean( dm )    
                    dmrescale=(dm/rescaledist)*100},
                   
                   {2  ## Create a distance matrix that is re-scaled
                    ## to the max distance
                    
                    rescaledist = max(dm);
                    dmrescale = (dm/rescaledist)*100}
                   )
        } else { dmrescale = dm }
    
    ## take the sd error as a unit for the measure
    combo = c(ts1,ts2)
    sdun = sd(dmrescale)
    mnun = median(dmrescale) ## the distance that gives us RR 50%
    
    ## sequence of radius from max to 0 in steps
    ## relative to SD
    radi = seq(mnun, 0, -(sdun/radiusspan))
    
    ## we take only the lower half, where it is certain that
    ## radius will produce RR < 25% 
    radi = radi[(length(radi)/2):length(radi)]
    
    
    kpt = ceiling(length(radi)/radiussample)
    rsamples = sample(1:kpt, 1)
    
    syssamp = seq(rsamples, rsamples + kpt*(radiussample-1), kpt)
    syssamp = syssamp[syssamp <= length(radi)]
  
    radi = radi[syssamp]
    
    ## delay and embed dimension to pass
    delay = del; embed = embdim;
    delay = 1; embed = 1

    optrad = vector()
                                        #  count = 0
    for (r in radi){
                                        #      count = count + 1
        ##  print( paste( "Checking radius:", r) )
        radius = r  
        res = crqa(ts1, ts2, delay, embed,
            rescale, radius, normalize, mindiagline,
            minvertline, tw, whiteline, recpt)
        
    #  print(res$rec)
        
        if (res$rec >= 2 & res$rec <= 5){
            optrad = r
            break}
        
    }
    
    if (length(optrad) == 0){
      optrad = NA}
    
    if(is.na(optrad)){
        warning("Optimal radius not found: try again choosing a wider radius span and larger sample size")}
    
    return ( list(radius = optrad, emddim = embdim, delay = del) )
    
}
