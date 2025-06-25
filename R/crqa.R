## originally written in R by Moreno I. Coco, 2013, (moreno.cocoi@gmail.com)
## crqa, inspired and adapted from a Matlab code developed at
## summer school of: Nonlinear Methods for Psychological Science
## organized by the University of Cincinnati, 2012
## most important update by Moreno I. Coco, 02/2019 using R code 
## by Sebastian Wallot (mdcrqa) v1.0, 13. April 2018

## arguments to pass to crqa:

## ts1, ts2: times series of integers indicating the states
## delay = nr. of lags
## embed = the embedding dimension, i.e., the lag intervals
## rescale = rescale the distance matrix before looking at the radius; 
##           if 1 (Mean Distance); if 2 (Max Distance), if 3 (Min Distance), if 4 (Euc Distance)
## radius =  distance to accept two points as recurrent (set it very 
##           small, if the series are categorical in nature)
## normalize = rescale the input variables for source data; 
##           if 1 (Unit interval); if 2 (z-score) 
## mindiagline = set a minimum diagonal line length
## mindiagline = set a minimum vertical line length
##  whiteline = FALSE # - flag to compute or not white vertical lines
##                    in the recurrence plot. Note, white lines are not
##                    yet used to derive any particular measure
##  recpt = FALSE # - flag to indicate whether the input ts1 is already
##                a recurrence plot

## tw = the size of the Theiler Window, the default is 0
## side = a string indicating whether the recurrence measures
## should be calculated in the "upper" triangle of the matrix
## "lower" triangle of the matrix, on the "whole" matrix
## method = a string vector indicating the type of recurrence analysis
##          options are: "rqa", "crqa" and "mdcrqa".
## metric = the distance measure to apply, 
##          default euclidean but see help(cdist) for more options

## datatype = (continuous, categorical) - nature of input data
##          default is continuous

## try below
## examples of categorical data 
## -- vectors
## ts1 = c("cat", "friend", "frenzy", "dog", "mum", "door")
## ts2 = c("miss", "shop", "dog", "mum", "incomprensible", "friend")
## matrices
## ts1 = cbind(c("cat", "friend", "frenzy", "dog", "mum", "door"),
##            c("miss", "shop", "dog", "mum", "incomprensible", "friend"))

## ts2 = cbind(c("friend", "frenzy", "dog", "mum", "door", "man"),
##            c("miss", "shop", "dog", "friend", "idea", "love"))

## examples of continuous data 
## -- vectors
# ts1 = c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1)
# ts2 = c(1, 1, 2, 2, 0, 0, 1, 2, 2, 2)
## -- matrices
# ts1 = cbind(c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1),
#       c(0, 0, 2, 1, 0, 1, 1, 2, 2, 1))
# ts2 = cbind(c(1, 1, 2, 2, 0, 0, 1, 2, 2, 1),
# c(2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 0, 0))

## starting parameters
## delay = 1; embed = 1; rescale = 1; radius = 0.001;
## normalize = 0; mindiagline = 2; minvertline = 2;
## tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
## method = 'mdcrqa'; metric = 'euclidean';  datatype = "continuous"

# crqa(ts2, ts1, delay, embed, rescale, radius, normalize, 
# mindiagline, minvertline, tw, whiteline, recpt, side, method, 
# metric, datatype)

# require(rdist) ## to choose the distance matrix choosing a specific metric
# require(Matrix) ## to manipulate sparse matrices 

.packageName <- 'crqa'

crqa <- function(ts1, ts2, delay = 1, embed = 1, rescale = 0,
                 radius = 0.001, normalize = 0, mindiagline = 2, minvertline = 2,
                 tw = 0, whiteline = FALSE, recpt = FALSE, side = "both", 
                 method = "rqa", metric = "euclidean", datatype = "continuous"){
  
  # print(data.frame(delay, embed, radius, rescale, 
  #                 normalize, mindiagline, minvertline, tw, whiteline, 
  #                 recpt, side, method, metric, datatype))
  
  ## first, need to check that the input variables whether all parameters have value
  # check input variables
  ## check if the input is a recurrence plot 
  if (recpt == FALSE){
    
    # first check whether input variables exist
    if (exists("ts1")) ts1 = ts1 else stop("No data has been specified for ts1")
    if (exists("ts2")) ts2 = ts2 else stop("No data has been specified for ts2")
    
    ## check if the method inputted is valid
    chkmet = method%in%c("rqa", "crqa", "mdcrqa")
    if (chkmet == F) stop("The method you have used is not valid")
    
    if (method == "rqa" | method == "crqa"){ ## data for rqa and crqa should be inputted as vector 
      if (is.matrix(ts1)) stop("Your data must consist of a single column of data.")  
      if (is.matrix(ts2)) stop("Your data must consist of a single column of data.")      
      
      ts1 = as.vector(as.matrix(ts1)) ## make sure data is a vector
      ts2 = as.vector(as.matrix(ts2))
      
      ## make sure the data is in a continuous format
      if (is.character(ts1) | is.factor(ts1) | is.character(ts1) | is.factor(ts1)){ 
        warning("Your input data was provided either as character or factor, and recoded as numerical identifiers")  
        tsnorm = numerify(ts1, ts2)
        ts1 = tsnorm$nwts1
        ts2 = tsnorm$nwts2
      }
      
      ## make sure they have the same length otherwise refer user
      if(length(ts1) != length(ts2)) stop("The input vectors have different length")
      ## check that the length of the series is not shorter than the phase embed*delay
      if (length(ts1) < embed*delay){ stop("Phase-space (embed*delay) longer than ts1")}  
      if (length(ts2) < embed*delay){ stop("Phase-space (embed*delay) longer than ts2")}      
    }
    
    if (method == "mdcrqa"){
      if (nrow(ts1) != nrow(ts2)) stop("ts1 and ts2 do not have the same number rows") 
      if (ncol(ts1) != ncol(ts2)) stop("ts1 and ts2 do not have the same number columns") 
      if (nrow(ts1) < (embed-1)*delay) stop("Insufficient number of data points to embedd time-series")
      ## make sure the data is in a continuous format
      if (is.character(ts1) | is.factor(ts1) | is.character(ts2) | is.factor(ts2)){ 
        warning("Your input data was provided either as character or factor, and recoded as numerical identifiers")  
        
        tsnorm = numerify(ts1, ts2)
        ts1 = tsnorm$nwts1
        ts2 = tsnorm$nwts2
      }
    }
    
    ## provide default values to all parameters if missing (R needs each exists to be in a single)
    if(exists("embed"))         embed = embed else embed = 1
    if(exists("delay"))         delay = delay else delay = 1
    if(exists("rescale"))       rescale = rescale else rescale = 0
    if(exists("normalize"))     normalize = normalize else normalize = 0
    if(exists("radius"))        radius = radius else radius = 1
    if(exists("mindiagline") & mindiagline > 2) mindiagline = mindiagline else mindiagline <- 2
    if(exists("minvertline") & minvertline > 2) minvertline = minvertline else minvertline <- 2
    if(exists("tw"))           tw = tw else tw = 0
    if(exists("whiteline"))    whiteline = whiteline else whiteline = F
    if(exists("recpt"))        recpt    = recpt      else recpt = F
    if(exists("side"))         side     = side       else side = "both"
    if(exists("method"))       method   = method     else method = "crqa"
    if(exists("metric"))       metric   = metric     else metric = "euclidean"
    if(exists("datatype"))     datatype = datatype   else datatype = "continuous"
    
    ##rescale the input data if necessary
    if (normalize > 0){
      switch (normalize,
              {1
                ## unit-interval
                ts1 = (ts1 - min(ts1));
                ts1 = ts1 / max(ts1);
                ts2 = (ts2 - min(ts2));
                ts2 = ts2 / max(ts2);},
              
              {2                     
                ## z-score                    
                ts1 = scale(ts1)      
                ts2 = scale(ts2)
              }
      )
    }
    
    ## check whether input data needs to be embedded 
    ## and need different procedures whether the data is vector (rqa/crqa) 
    ## or a matrix (mdcrqa)
    
    if (embed > 1) {
      if (method == 'rqa' | method == 'crqa'){
        newLength <- length(ts1) - (embed-1)*delay
        tempTs1 <- ts1[1:newLength]
        for (i in seq(2,embed)) {
          tempTs1 <- cbind(tempTs1, ts1[(1+(delay*(i-1))):(newLength+delay*(i-1))])
        }
        ts1 <- tempTs1
        rm(tempTs1)
        tempTs2 <- ts2[1:newLength]
        for (i in seq(2,embed)) {
          tempTs2 <- cbind(tempTs2, ts2[(1+(delay*(i-1))):(newLength+delay*(i-1))])
        }
        ts2 <- tempTs2
        rm(tempTs2)
      }
      
      if (method == 'mdcrqa'){
        newLength <- dim(ts1)[1] - (embed-1)*delay
        tempTs1 <- ts1[1:newLength,]
        for (i in seq(2,embed)) {
          tempTs1 <- cbind(tempTs1,ts1[(1+(delay*(i-1))):(newLength+delay*(i-1)),])
        }
        ts1 <- tempTs1
        rm(tempTs1)
        tempTs2 <- ts2[1:newLength,]
        for (i in seq(2,embed)) {
          tempTs2 <- cbind(tempTs2,ts2[(1+(delay*(i-1))):(newLength+delay*(i-1)),])
        }
        ts2 <- tempTs2
        rm(tempTs2)
      }
    }
    
    ## just to have the length of matrix saved
    dm <- as.matrix(cdist(ts1, ts2, metric = metric))
    
    ## Find indeces of the distance matrix that fall
    ## within prescribed radius.
    if (rescale > 0){
      switch(rescale,
             {1  ## Create a distance matrix that is re-scaled
               ## to the mean distance
               
               rescaledist = mean(dm)    
               dmrescale   = dm/rescaledist},
             
             {2  ## Create a distance matrix that is re-scaled
               ## to the max distance
               
               rescaledist = max(dm);
               dmrescale   = dm/rescaledist},
             {3 ## Create a distance matrix that is rescaled 
               ## to the min distance
               rescaledist = min(dm);
               dmrescale   = dm/rescaledist},
             
             { 4 ## Create a distance matrix that is rescaled 
               ## to the euclidean distance
               dmrescale <- dm/abs(sum(dm)/(nrow(dm)^2-nrow(dm)))}
      )
    } else { dmrescale = dm }
    ## Compute recurrence matrix
    
    v1l = nrow(dmrescale); v2l = ncol(dmrescale) ## save the dimension of the matrix
    
    ind = which(dmrescale <= radius, arr.ind = TRUE);
    r = ind[,1]; c = ind[,2]
    
  } else { ## take as input an RP directly
    if (exists("ts1")) ts1 = ts1 else stop("No data has been specified for ts1")
    ## as usual R needs fiddly code to make sure about identify of data
    ts1 = matrix(as.logical(ts1), ncol = ncol(ts1))
    v1l = nrow(ts1); v2l = ncol(ts1)        
    
    ## matrix needs to be logical
    ind = which(ts1 > 0, arr.ind = TRUE)
    
    ## just a trick to reduce the number of lines
    ## of the code
    r = ind[,1]; c = ind[,2]
    
  }
  
  if (length(r) != 0 & length(c) != 0){ ##avoid cases with no recurrence
    
    S = sparseMatrix(r, c, dims = c(v1l, v2l))
    ## this is the recurrent plot
    ## transpose it to make identical to Marwan
    S = t(S)
    
    ## apply the theiler argument here to recurrence matrix
    ## Marwan blanks out the recurrence along the diag
    S = theiler(S, tw)
    
    if (side == "upper"){
      ## if only the upper side is of interest
      ## it blanks out the lowest part
      S = as.matrix(S)
      S[lower.tri(S, diag = TRUE)] = 0
      S = Matrix(S, sparse = TRUE)
    }
    
    if (side == "lower"){
      ## viceversa
      S = as.matrix(S)
      S[upper.tri(S, diag = TRUE)] = 0
      S = Matrix(S, sparse = TRUE)
    }
    
    if (side == "both"){
      ## just keep it as is.
      S = S}
    
    spdiagonalize = spdiags(S) ##  spdiags should have decent speed 
    B = spdiagonalize$B
    
    ##calculate percentage recurrence by taking all non-zeros
    
    numrecurs = length(which(B == TRUE));
    percentrecurs = (numrecurs/((v1l*v2l)))*100;
    
    ####################################################################
    ####################################################################
    
    ## Computing the line counts
    
    ## This section finds the index of the zeros in the matrix B,
    ## which contains the diagonals of one triangle of the
    ## recurrence matrix (the identity line excluded).
    
    ## The find command indexes the matrix sequentially
    ## from 1 to the total number of elements.
    ## The element numbers for a 2X2 matrix would be [1 3; 2 4].
    ## You get a hit for every zero. If you take the difference
    ## of the resulting vector, minus 1, it yields the length of an
    ## interceding vector of ones, a line. Here is an e.g.
    ## using a row vector rather than a col. vector, since it types
    ## easier: B=[0 1 1 1 0], a line of length 3.
    ## find( B == 0 ) yields [1 5], diff( [1 5] ) -1 = 3,
    ## the line length.
    ## So this solution finds line lengths in the interior of
    ## the B matrix, BUT fails if a line butts up against either
    ## edge of the B matrix, e.g. say  B = [0 1 1 1 1],
    ## which( B == 0) returns a 1, and you miss the line of length 4.
    ## A solution is to "bracket" B with a row of zeros at each
    ## top and bottom.
    
    ## Bracket B with zeros
    if (is.vector(B)) {
      false = rep(FALSE, length(B)) ##cases where B is a vector
      B = rbind(false, B, false, deparse.level = 0)
    } else  {
      false = rep(FALSE, ncol(B))
      B = as.matrix(B)
      ## need to transform the sparseMat into normal to bracket it
      B = rbind(false, B, false, deparse.level = 0)
    }
    
    ## Get list of line lengths, sorted from largest to smallest
    diaglines = sort( diff(which(B == FALSE) ) -1, decreasing = TRUE)
    
    ## Delete line counts less than the minimum diagonal.
    diaglines = diaglines[diaglines >= mindiagline]
    ## diaglines(diaglines>200)=[]; # Can define a maximum line length too.
    
    ## exlude the rare cases where there are no diaglines
    
    if(length(diaglines) != 0){
      
      numdiaglines = length(diaglines) ## extract the length of diag
      maxline = max(diaglines)
      meanline = mean(diaglines)
      
      tabled = as.data.frame(table(diaglines))
      
      total = sum(tabled$Freq)       
      p = tabled$Freq/total
      
      ##remove zero probability..it should not be necessary
      del = which(p == 0 )
      if (length(del) > 0) {
        p = p[-del]
      }
      
      ## entropy log2, and relative entropy divided by max
      entropy = - sum(p*log(p))    
      relEntropy = entropy/(-1*log(1/nrow(tabled)))
      
      ## entropy/max entropy: comparable across contexts and conditions.
      
      pdeter = sum(diaglines)/numrecurs*100
      ## percent determinism: the predictability of the dynamical system 
      
      ## calculate laminarity and trapping time
      restt = tt(S, minvertline, whiteline)            
      lam = restt$lam; TT = restt$TT; max_vertlength = restt$max_vertlength
      
      ## let's calculate categorical entropy
      if (side == 'both' & datatype == 'categorical' & radius <= .1){ 
        ## we need a full RP and data has to be categorical
        ## we need to input directly the indeces of the 
        ## recurrence plot and the size of the matrix
        size   = dim(S)
        catH   = catEnt(ind, size)
      } else {
        catH = NA
      }
      
    } else {
      
      numdiaglines = 0; maxline = 0; pdeter = NA;
      entropy = NA; relEntropy = NA; meanline = 0
      lam = 0; TT = 0; catH = NA; max_vertlength = NA; 
      RP = NA; 
    }
    
    results = list(RR = percentrecurs, DET = pdeter, 
                   NRLINE = numdiaglines, maxL = maxline, 
                   L = meanline, ENTR = entropy, 
                   rENTR = relEntropy,
                   LAM = lam, TT = TT, catH = catH, 
                   max_vertlength = max_vertlength, RP = S)
    
  } else { # print (paste ("No recurrence found") )
    results = list(RR = 0, DET = NA, NRLINE = 0,
                   maxL = 0, L = 0, ENTR = NA, rENTR = NA,
                   LAM = NA, TT = NA, catH = NA, 
                   max_vertlength = NA, RP = NA)}  
  
  return (results)
  
}




