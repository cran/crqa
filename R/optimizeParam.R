## GNU > 2 License:
## written by Moreno I. Coco (moreno.cocoi@gmail.com) and James Dixon

## Iterative procedure exploring a combination of parameter values
## to obtain maximal recurrence between two time-series
## The procedure should be applied on several trials to find the
## "optimal" setting for the parameters...
## it should be applied on time-series with continuous data

## TODO: probably, for windows the selection of size/embedding
## has to be done at the same time. 

## some tips

## Recurrence measures change with changing parameters:
## Explore parameters for both signals same way as for
## phase space reconstruction / auto recurrence.
## If you find the same parameters for both signals: Great!
## If you find different parameters for each signal
## (which is more common):
## Run Cross RQA once for each set of parameters
## and make sure the results will pattern the same way
## across both parameter settings
## If so, you can be more confident that any results are
## not due to parameter differences.

## maxlag = 20 ## maximum lag to check for average mutual information

## radiusspan = 2 # increasing radius step means smaller steps
# #radiussample = 10 ## number of radius points within the steps
## to be sampled

## min.rec = minimal recurrence value
## max.rec = maximal recurrence value

#ts1 = runif(100)
#ts2 = runif(100)

.packageName <- 'crqa'

optimizeParam <- function(ts1, ts2, par, min.rec = 2, max.rec = 5){
  
  ## Initialize attached variables
  method = metric = maxlag = radiusspan = 
  normalize = rescale =  mindiagline = minvertline = tw =
  whiteline = recpt = side = datatype =  fnnpercent = 
  typeami =  nbins = criterion = threshold =
  maxEmb = numSamples = Rtol = Atol =  NULL
  
  ## Assign parameters
  for (v in 1:length(par)) assign(names(par)[v], par[[v]])
  
  
  if (method == 'rqa' | method == 'crqa'){
    ## we estimate the paramaters for radius and embedding considering 
    ## two one-dimensional time-series
    
    if (is.vector(ts1) != TRUE){ ts1 = as.numeric(as.matrix(ts1) ) }
    if (is.vector(ts2) != TRUE){ ts2 = as.numeric(as.matrix(ts2) ) }
    
    ##  A: Choose a delay that will accommodate both ts.
    ##  For example if one has a considerably longer delay
    ##  indicated than another, you may need to select the longer
    ##  delay of the two because that ensures that new information
    ##  is gained for both by using the higher delay.
    ##  If delays are close to each other, you may want to use a
    ##  delay somewhere in between the two
    
    mi1 = as.numeric(mutual(ts1, lag.max = maxlag, plot = FALSE))
    mi2 = as.numeric(mutual(ts2, lag.max = maxlag, plot = FALSE))
    ##  average mutual information    
    mi = ami(ts1, ts2, 1:maxlag)
    
    ## global minimum for each series and ami, given lag max
    m1 = min(mi1); m2 = min(mi2); m12 = min(mi)
    mis = c(m1, m2, m12)
    
    ## find the minimum dip across all (mi1, mi2, ami(1/2))
    
    if (typeami == "mindip"){
      minmi = which(mis == min(mis))
      
      if (minmi == 1)  lag = which(mi1 == m1)
      if (minmi == 2)  lag = which(mi2 == m2)
      if (minmi == 3)  lag = which(mi == m12)
      
      maxlag = lag
    }
    
    ## just double check that the maximum lag is the same
    ## with the minimum dip, if not just take the longest lag
    
    if (typeami == "maxlag"){
      
      lg1 = which(mi1 == m1)
      lg2 = which(mi2 == m2)        
      lg12 = which(mi == m12)
      
      lags = c(lg1, lg2, lg12)
      
      maxlag = lags[which(lags == max(lags))]
    }
    
    if (length(maxlag) > 1){
      ## just to address those cases more than a lag has the same max/min
      maxlag = maxlag[sample(1:length(maxlag),1)]
    }
    
    del = maxlag - 1 
    
    ## another solution would be to use the global minimum
    ## but it could result into longer delays being selected
    
    ## B: Determine embedding dimension: FNN function bottoms out
    ## (it buys nothing to add more dimensions).
    ## If the embedding dimension indicated for the two files is
    ## different, you want to go with the higher embedding
    ## dimension of the two to make sure that you have
    ## sufficiently unfolded both time series.
    ## It is generally safe to overestimate embedding dimensions.
    ## Though this can amplify the influence of noise for
    ## high embedding dimensions. 
    
    ## TODO: bring out this parameters to optimize them also    
    ## m	 maximum embedding dimension = 20
    ## d	 delay parameter = estimate at 1
    ## t	 Theiler window = 0
    ## in cross-recurrence LOC is important
    ## rt	 escape factor = leave it default
    ## eps	 neighborhood diameter = leave default
    
    embdts1 = false.nearest(ts1, m = maxEmb, d = del, t = 0, rt = 10,
                            eps = sd(ts1)/10)
    
    ## get a percentage of reduction of false neighbours
    ## based on the first dimension
    
    fnnfraction1 = embdts1[1,]
    fnnfraction1 = fnnfraction1[which(is.na(fnnfraction1) == FALSE)]
    emdthd1 = fnnfraction1[1]/fnnpercent
    
    emdix1 = which(diff(fnnfraction1) < -emdthd1)
    
    if (length(emdix1) == 1){
      emdmints1 = as.numeric(emdix1) + 1
    } else if (length(emdix1) > 1){
      ## there is no gain when embedding
      emdmints1 =  as.numeric(tail(emdix1,1) + 1)
    }  else { 
      emdmints1 = 1
    }
    
    ## to adjust for the diff 
    
    ## alternative method, just get the minimum
    ## at the moment the method is commented.
    # embmints1 = as.numeric( which(fnnfraction1 == min(fnnfraction1)))
    
    embdts2 = false.nearest(ts2, m = maxEmb, d = del, t = 0, rt=10,
                            eps=sd(ts2)/10)
    
    fnnfraction2 = embdts2[1,]
    fnnfraction2 = fnnfraction2[which(is.na(fnnfraction2) == FALSE)]
    emdthd2 = fnnfraction2[1]/fnnpercent
    
    emdix2 = which(diff(fnnfraction2) < -emdthd2)
    
    if (length(emdix2) == 1){
      emdmints2 = as.numeric(emdix2) + 1
    } else if (length(emdix2) > 1){
      emdmints2 =  as.numeric(tail(emdix2,1) + 1)
    }  else {         ## there is no gain when embedding
      emdmints2 = 1
    }
    
    # embmints2 = as.numeric( which(fnnfraction2 == min(fnnfraction2)))
    
    if ( length(emdmints1) > 1){ emdmints1 = emdmints1[1] }
    if ( length(emdmints2) > 1){ emdmints2 = emdmints2[1] }
    
    embdim = max( c(emdmints1, emdmints2) )        
    
    ## the radius could be theoretically estimated also for 
    ## the multidimensional case, 
    ## but only if the variables are z-scored or are coming from the same measuremnt 
    
    ## C: Use the optimal parameters of embedding and delay
    ## from above and then explore an optimal radius
    
    ## first we need to see in the rescaled matrix what is
    ## the range of values, so that the radius interval
    ## can be consciously estimated.
  }
  
  if (method == 'mdcrqa'){
    ## we use the functions for multidimensional estimation
    ## do it separately for ts1 and ts2, then get the largest 
    ## value between the two
    
    
    ## estimate delay
    tryCatch( { delay_ts1 = mdDelay(ts1, nbins, maxlag, criterion, threshold)}
              , error = function(e) {    
                stop("Delay not found, try increase the threshold, 
or explore more lags(maxlag) on a larger number of bins (nbins)")})
    
    tryCatch( { delay_ts2 = mdDelay(ts2, nbins, maxlag, criterion, threshold)}
              , error = function(e) {    
                stop("Delay not found, try increase the threshold, 
or explore more lags (maxlag) on a larger number of bins (nbins)")})
    
    del = max(c(delay_ts1, delay_ts2))
    
    ## get a percentage of reduction of false neighbours
    ## based on the first dimension
    fnn_ts1   = mdFnn(ts1, delay_ts1, maxEmb, numSamples, Rtol, Atol)  
    fnn_ts2   = mdFnn(ts2, delay_ts2, maxEmb, numSamples, Rtol, Atol)
    
    fnnfraction1 = fnn_ts1[[1]]
    ixi_min = which(fnnfraction1 == min(fnnfraction1))
    if (length(ixi_min) > 1){
      embmints1 = ixi_min[1]
    } else {
      embmints1 = ixi_min
    }
    
    
    fnnfraction2 = fnn_ts2[[1]]
    ixi_min = which(fnnfraction2 == min(fnnfraction2))
    if (length(ixi_min) > 1){
      embmints2 = ixi_min[1]
    } else { 
      embmints2 = ixi_min
      }
    
    embdim = max( c(embmints1, embmints1) ) 
    
  }
  
  if (radiusspan <= 1){
    stop("Radius span too small, please choose a larger value")
  }
  
  
  dm <- as.matrix(cdist(ts1, ts2, metric = metric))
  ## Find indeces of the distance matrix that fall
  ## within prescribed radius.
  if (rescale > 0){
    switch(rescale,
           {1  ## Create a distance matrix that is re-scaled
             ## to the mean distance
             
             rescaledist = mean(dm, na.rm = T)    
             dmrescale   = dm/rescaledist},
           
           {2  ## Create a distance matrix that is re-scaled
             ## to the max distance
             
             rescaledist = max(dm, na.rm = T);
             dmrescale   = dm/rescaledist},
           {3 ## Create a distance matrix that is rescaled 
             ## to the min distance
             rescaledist = min(dm, na.rm = T);
             dmrescale   = dm/rescaledist},
           
           { 4 ## Create a distance matrix that is rescaled 
             ## to the euclidean distance
             dmrescale <- dm/abs(sum(dm)/(nrow(dm)^2-nrow(dm)))}
    )
  } else { dmrescale = dm }
  
  
  sdun = sd(dmrescale)
  mnun = median(dmrescale)*3
  ## Multiplier to increase range of candidate radii 
  
  ## generate a sequence of radius
  radi = seq(mnun, 0, -(sdun/radiusspan))
  ##    radi = radi[(length(radi)/2):length(radi)]
  ##    kpt = ceiling(length(radi)/radiussample)
  ##    rsamples = sample(1:kpt, 1)
  
  ##    syssamp = seq(rsamples, rsamples + kpt * (radiussample - 1), kpt)
  ##    syssamp = syssamp[syssamp <= length(radi)]
  ##    radi = radi[syssamp]
  delay = del
  embed = embdim
  optrad = vector()
  
  end.flag <- 0
  
  while(end.flag == 0){
    ## location with largest radius
    hi.loc <- 1
    ## location with smallest radius
    lo.loc <- length(radi)     
    
    #print(c("hi.loc",hi.loc))
    #print(c("lo.loc",lo.loc))
    ## take middle location as the current radius to be tested 
    
    curr.loc <- round(length(radi)/2)
    radi[curr.loc] -> r
    # print(radi)
    
    #print(c("r",r))
    #print(c("curr.loc",curr.loc))
    #print(c("embed",embed))
    #print(c("delay",delay))
    
    res = crqa(ts1, ts2, delay, embed, rescale, r, normalize, 
               mindiagline, minvertline, tw, whiteline, recpt, side, 
               method, metric, datatype) 
    # print(c("recurr",res$RR))
    ## print("#######")
    ## print("#######")       
    
    if (res$RR >= min.rec & res$RR <= max.rec) {
      optrad = r
      ## print(c("radius",optrad))
      ## print(c("embed",embed))
      ## print(c("delay",delay))
      end.flag<-1
    } else {
      
      if (res$RR < min.rec){
        ## if rec less than min rec, curr.loc becomes new lo	
        curr.loc -> lo.loc
      }
      if (res$RR > max.rec){
        ## if rec greater than max rec, curr.loc becomes new hi	
        curr.loc -> hi.loc
      }   	 
      if((lo.loc-hi.loc) < 2){
        ## if less than 2 radi locs remaining, no optimal radius	
        end.flag<-1
        warning("Optimal Radius Not found: try again choosing a wider radius span and larger sample size")
      }
    }
    
    ## replace radi vector with remaining unsearched vector
    radi <- radi[hi.loc:lo.loc]
  } ## end while         
  
  if (length(optrad) == 0) {
    optrad = NA
  }    
  
  if (!is.na(optrad)){
    return(list(radius = optrad, emddim = embdim, delay = del))
  }
  
}
