#MDFNN The mdfnn function computes the percentage of false nearest
#   neighbors for multidimensional time series as a function of embedding
#   dimension based on Kennel, M. B., Brown, R., & Abarbanel, H. D. (1992).
#   Determining embedding dimension for phase-space reconstruction using a
#   geometrical construction. Physical review A, 45, 3403.
#
#   Inputs:
#   Required arguments:
#    data = an m*n matrix, were m is the number of data points and n is the
#    number of dimensions of the time series.
#
#    tau = time delay for embedding.
#
#   Optional arguments:
#    maxEmb = maximum number of embedding dimensions that are considered.
#     Default is 10.
#
#    doPlot = plot false nearest neighbors if 1, supress plot if 0.
#     Default is 1.
#
#    noSamples = number of randomly drawn coordinates from phase-space used
#     to estimate the percentage of false-nearest neighors.
#     eflaut is 500.
#
#    Rtol = First distance criterion for separating false neighbors.
#     Default is 10.
#
#    Atol = Second distance criterion for separating false neighbors.
#     Default is 2.
#
#   Outputs:
#    fnnPerc = percentage of false neighbors for each embedding.
#  
#    embTimes = number of times the multidimensional time series was
#    embedded with delay tau. Note, that embTimes = 1 means no embedding,
#    and embTimes = 2 means embedding the multidimensional time series
#    once. Hence, the resulting falase neibor percentages result from
#    embeddeding data (embTimes - 1) times
#
#   Matlab Version: 1.0, 22 June 2018
#   Authors:
#     Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
#     Dan Moenster, Aarhus University
#
#   Reference:
#     Wallot, S., \& M{\o}nster, D. (2018). Calculation of average mutual
#     information (AMI) and false-nearest neighbors (FNN) for the
#     estimation of embedding parameters of multidimensional time-series in
#     Matlab. Front. Psychol. - Quantitative Psychology and Measurement
#     (under review)

# recoded in R by Moreno I. Coco (moreno.cocoi@gmail.com), 25/02/2019

mdFnn = function(data, tau = 1, maxEmb = 10, numSamples = 500, Rtol = 10, Atol = 2){
  
  ## check the data type and possible errors in it
  tdata = class(data)[1]
  
  if (tdata == "data.frame"){data = as.matrix(data)} ## convert data.frames into matrices
  
  if (sapply(data, is.numeric)[1] != TRUE){
    stop('Input is not numeric')}
  
  if (is.vector(data) == TRUE & length(data) <= 1){
    stop('Input must be a vector or matrix')}
  
  if (is.matrix(data) == TRUE & ncol(data) <= 1){
    stop('Input must be a vector or matrix')}
  
  
  # Now do the actual work to find FNN
  #
  maxEmb = maxEmb+1
  
  if (is.vector(data) == TRUE){ dims = 1; N = length(data) # get dimensionality of time series
  } else {
    dims = ncol(data); N = nrow(data) 
  }
  
  fnnPerc = rep(0,(maxEmb-1)) #fnnPerc = 100       # first FNN
  
  Ra = sum(diag(var(data))) # estimate of attractor size
  
  if ((N - tau *(maxEmb-1)) < numSamples){ # check whether enough data points exist for random sampling
    numSamples = N - tau*(maxEmb-1)
    samps = 1:numSamples
  } else {
    samps = 1:(N - tau*(maxEmb-1)) # generate random sample
    samps = samps[sample(length(samps), numSamples)]
  }
  
  embData = vector()
  #i = 2
  for (i in 1:maxEmb){ # embed data
    # print(i)
    embData = cbind(embData, data[(1 + (i-1)*tau):(N-(maxEmb-i)*tau), 1:dims], deparse.level = 0)
    dim(embData)
    embData[1:10, ]
    dists = pdist2(embData^2, embData^2);
    dim(dists)
    dists[1:10, 1:10]
    r2d1 = yRd1 = 0
    
    j = 1
    for (j in 1:numSamples){ # get nearest neighbors and distances
      m = sort(dists[, samps[j]], index.return = TRUE);
      
      temp  = m$x
      coord = m$ix
      
      r2d1[j] = temp[2]
      yRd1[j] = coord[2]
    }
    
    r2d1[1:10]
    yRd1[1:10]
    
    #r2d = r2d1 #!
    #yRd = r2d1 #!
    
    if (i == 1){ #!
      r2d = r2d1
      yRd = r2d1
    } else {
      
      fnnTemp = 0
      j = 1
      for (j in 1:length(r2d1)){
        temp = dists[, samps[j]] 
        length(temp)
        temp[1:10]
        
        # check whether neighbors are false
        dst = sqrt((temp[yRd1[j]] - r2d[j])/r2d[j]) ## need to store it outside because it could return NAs
        
        if (is.na(dst) == F){ ## need to check whether dst was created
          
          if (dst  > Rtol || abs(temp[yRd1[j]] - r2d[j])/Ra > Atol){
            
            fnnTemp[j] = 1;
          } else {
            fnnTemp[j] = 0;
          }
        }
      }
      
      fnnPerc[i-1] = 100*sum(fnnTemp)/length(fnnTemp); # compute percentage of FNN
      r2d = r2d1;
      yRd = r2d1;
      
    }
  }
  
  embTimes = 1:maxEmb;
  return(list(fnnPerc = fnnPerc, embTimes = embTimes))
}