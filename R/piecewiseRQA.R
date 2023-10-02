## piece-wise RQA
# typeRQA = c("full", "diagonal") ## TODO for windowed method

piecewiseRQA <- function(ts1, ts2, blockSize, delay = 1, embed = 1, rescale = 0,
                         radius = 0.001, normalize = 0, mindiagline = 2, minvertline = 2,
                         tw = 0, whiteline = FALSE, recpt = FALSE, 
                         side = "both", method = "crqa", metric = "euclidean", 
                         datatype = "continuous", typeRQA = "full", 
                         windowsize = NA){
  
  if (exists("ts1")) ts1 = ts1 else stop("No data has been specified for ts1")
  if (exists("ts2")) ts2 = ts2 else stop("No data has been specified for ts2")
  
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
  
  ## initialise all parameters value if missing (R needs this weird repetition to be the code in one line)
  if(exists("embed"))         embed = embed else embed = 1
  if(exists("delay"))         delay = delay else delay = 1
  if(exists("normalize"))     normalize = normalize else normalize = 0
  if(exists("radius"))        radius = radius else radius = 1
  if(exists("mindiagline") & mindiagline > 2) mindiagline = mindiagline else mindiagline <- 2
  if(exists("minvertline") & minvertline > 2) minvertline = minvertline else minvertline <- 2
  if(exists("typeRQA"))     typeRQA = typeRQA else typeRQA = "full"
  
  ## i guess these parameters should be assigned depending on the ts1, and ts2
  nrows = max(1:length(ts1) - delay*(embed-1))
  ncols = max(1:length(ts2) - delay*(embed-1))
  
  # create the big RP to fill in like a mosaic with smaller RPs 
  
  piecewiseRP = Matrix(0, nrow = nrows, ncol = ncols, sparse = TRUE);
  
  nrs = seq(blockSize, length(ts1) - delay * (embed-1), blockSize)
  nrc = seq(blockSize, length(ts2) - delay * (embed-1), blockSize)
  
  # loop through time-series in window sizes and create sub-RPs and sub-CRPs
  i = nrs[1]; j = nrc[2]
  for (i in nrs) {
    for (j in nrc){
      
      win1 = ts1[(1 + (i - blockSize)):(i + delay * (embed - 1))]
      win2 = ts2[(1 + (j-blockSize)):(j + delay * (embed-1))]
      
      temp =  crqa(win1, win2, delay, embed, rescale,
                   radius, normalize, mindiagline, minvertline,
                   tw, whiteline, recpt, side, method, metric,
                   datatype)
      temp = temp$RP
      
      piecewiseRP[(1 + (j - blockSize)):j, (1 + (i-blockSize)) : i] = temp;
      
    }  
  }
  
  
  # if the number of (embedded) data points is not a multiple of the window size, 
  # calculate sub-RPs and sub-CRPs for the edges of the piecewise RP
  
  ## compute the edges
  er = length(ts1)-delay*(embed-1)
  ec = length(ts2)-delay*(embed-1)
  
  if (length(ts1) - i - delay*(embed-1) > 0 | length(ts2) - j - delay*(embed-1) > 0){ 
    
    for (i in er){ # is this always just a single number?
      for (j in nrc){
        
        win1 = ts1[(1 + (i - blockSize)):(i + delay * (embed - 1))]
        win2 = ts2[(1 + (j-blockSize)):(j + delay * (embed-1))]
        
        temp =  crqa(win1, win2, delay, embed, rescale,
                     radius, normalize, mindiagline, minvertline,
                     tw, whiteline, recpt, side, method, metric,
                     datatype)
        temp = temp$RP
        
        piecewiseRP[(1 + (j - blockSize)):j, (1 + (i-blockSize)) : i] = temp;
        
      }
      
    }
    
    for (i in nrs){
      for (j in ec){
        
        win1 = ts1[(1 + (i - blockSize)):(i + delay * (embed - 1))]
        win2 = ts2[(1 + (j-blockSize)):(j + delay * (embed-1))]
        
        temp =  crqa(win1, win2, delay, embed, rescale,
                     radius, normalize, mindiagline, minvertline,
                     tw, whiteline, recpt, side, method, metric,
                     datatype)
        
        temp = temp$RP
        
        piecewiseRP[(1 + (j - blockSize)):j, (1 + (i-blockSize)) : i] = temp;
        
        
      }
    }
    
    win1 = ts1[(1 + (er - blockSize)):(er + delay * (embed - 1))]
    win2 = ts2[(1 + (ec - blockSize)):(ec + delay * (embed-1))]
    
    temp =  crqa(win1, win2, delay, embed, rescale,
                 radius, normalize, mindiagline, minvertline,
                 tw, whiteline, recpt, side, method, metric,
                 datatype)
    
    temp = temp$RP
    
    piecewiseRP[(1 + (ec - blockSize)):ec, (1 + (er - blockSize)) : er] = temp;
    
  }
  
  ## once we have computed the piecewise RP we need to extract the measures
  ## that we with to extract in it
  
  if (typeRQA == "full"){ ## we want to get the classic RQA measures
    recpt = TRUE
    
    results =  crqa(piecewiseRP, ts2 = NA, delay, embed, rescale,
                    radius, normalize, mindiagline, minvertline,
                    tw, whiteline, recpt, side, method, metric,
                    datatype)
  }
  
  
  if (typeRQA == "diagonal"){ ## we want to get the classic RQA measures
    
    RP = piecewiseRP
    
    if (is.logical(RP) != T){ # we have some point that recur
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
    
    results = list(profile = drpd, maxrec = maxrec, maxlag = maxlag)
    
  }

  return(results)
  
}


