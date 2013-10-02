# GNU > 2 License:
## written by Moreno I. Coco (moreno.cocoi@gmail.com)

# Iterative procedure exploring a combination of parameter values
# to obtain maximal recurrence between two time-series
# The procedure should be applied on several trials to find the
# "optimal" setting for the parameters...
# it should be applied on time-series with continuous data

# TODO: probably, for windows the selection of size/embedding
# has to be done at the same time. 


# source("./ami.R")       ## average mutual information
## between 2 series
# source("./crqa.R")

# ts1 = read.table("./datatype/ContSP1.txt", header = FALSE)
# ts2 = read.table("./datatype/ContSP2.txt", header = FALSE)
# lgM = 100 ## maximum lag to check for average mutual information
# steps = seq(1, 10, 1) ##how many points we should look ahead to define a local minimum
  
# cut.del = seq(1, 40,1)  ## cut off of delay is set at 500 ms (25ms * 20 delay)

.packageName <- 'crqa'

optimizeParam <- function(ts1, ts2, par){
  rescale = normalize = mindiagline = minvertline =
      lgM = steps = cut.del = recpt =  whiteline = NULL
  ## initialize attached variables

  for (v in 1:length(par)) assign(names(par)[v], par[[v]]) ## assign parameters
  
  #  require(tseriesChaos)

  ## check that the input is vector and numeric
  
  if (is.vector(ts1) != TRUE){ ts1 = as.numeric(as.matrix(ts1) ) }
  if (is.vector(ts2) != TRUE){ ts2 = as.numeric(as.matrix(ts2) ) }


  mi1 = as.numeric( mutual(ts1, lag.max = lgM, plot = FALSE))
  mi2 = as.numeric( mutual(ts2, lag.max = lgM, plot = FALSE))
        
 ## iterate over a range of possible steps to find out the optimal one
        
  lag1 = lag2 = vector()
  for (s in steps){
    
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
        
## for each time-series take the lag where MI is minimal
  
  ix.lg1 = which( mi1[lag1] == min(mi1[lag1]) )
  if (length(ix.lg1) > 1){
    ix.lg1 = sample(ix.lg1, 1)
    lg1 = lag1[ix.lg1]}

  ix.lg2 = which( mi2[lag1] == min(mi2[lag2]) )
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
#    print(delay)
          
    if (lg1 > lg2 & abs(lg1 - lg2) > delay){
#    print( c(delay, lg1) )  
      delay.in = c(delay.in, lg1)
    }
    
    else if (lg2 > lg1 & abs(lg2 - lg1) > delay){
#    print( c(delay, lg2) )            
      delay.in = c(delay.in, lg2)
    }
    
    else {
      delay.in = c(delay.in, mean(c(lg1, lg2)))
    }
          
#   print(delay.in)
    
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
# t	 Theiler window = 0 ## it seems a sensible thing to do
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

  ## take the sd error as a unit for the measure
  combo = c(ts1,ts2)
  sdun = sd(combo)
  mnun = min(combo[which(combo != 0)])
  ## mxun = max(combo)

  ## TODO: make this variable an argument of the function
  granradi = 15 ## times smaller of the SD unit 
  ## to create the observational step
        
  ## we start with a max radius larger as
  ## one unit of SD
        
  radi = seq(mnun, 0, -(sdun/granradi)) 

  ## initialize arguments for CRQA
  ## parameters for the window analysis
  delay = del; embed = embdim;
  rescale = 1; normalize = 2; minline = 2
  
  for (r in radi){
      ##  print( paste( "Checking radius:", r) )
      radius = r  
      res = crqa(ts1, ts2, delay, embed,
          rescale, radius, normalize, mindiagline,
          minvertline, whiteline, recpt)
      
      if (res$rec < 5){
          break}
      
  }
  
  return ( list(radius = r, emddim = embdim, delay = del) )
  
}












