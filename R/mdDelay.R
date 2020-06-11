# MDDELAY Estimates time delay for embedding of multivariate times series.
#    The function plots the mutual information for multivariate times series
#   data, so the user can estimate the optimal value of the time delay for
#   embedding the data. The function also returns an estimate of the
#   optimal time delay, using simple methods, such as the mean of the lag
#   for which the auto mutual information for each of the variables
#   (columns) is less than a threshold, such as 1/e.
#
#   This function currently implements the uniform multivariate embedding
#   method.

#   Inputs:
#   Required arguments:
#     data - a matrix with multiple timeseries, one in each column.
#
#   Optional arguments:
#     maxlag: The maximum time lag for which AMI is computed. Default = 10.
#
#     nbins: The number of bins used to construct the histograms for
#       computing AMI. Default = 10.
#
#     criterion: The criterion used for finding the optimal delay. Possible
#     values are:
#       'firstBelow' to use the lowest delay at which the AMI
#         function drops below the value set by the threshold parameter.
#       'localMin' to use the position of the first local minimum of the
#         AMI function.
#       Default: 'firstBelow'
#
#     threshold: The threshold value to select the delay when AMI drops
#       below threshold. Default = exp(-1)
#
#     plottype: Determines how the AMI is plotted. Possible values are
#     'mean', 'all', 'both', 'none'. Default = 'mean'
#
#   Outputs:
#     delay: The estimated optimal delay given the input and critera.
#
#   Based on Matlab Version: 1.0, 22 June 2018
#   Authors:
#   Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
#     Dan Moenster, Aarhus University
#
#   Reference:
#     Wallot, S., \& M{\o}nster, D. (2018). Calculation of average mutual
#     information (AMI) and false-nearest neighbors (FNN) for the
#     estimation of embedding parameters of multidimensional time-series in
#   Matlab. Front. Psychol. - Quantitative Psychology and Measurement
#     (under review)

# R translation by Moreno I. Coco (moreno.cocoi@gmail.com)
## arguments are: nbins = 10 Default
##                maxlag  = 10 Default 
##                criterion = 'firstBelow' 'localMin'
##                threshold = exp(-1) Default
## nbins = 10; maxlag = 10; criterion = "firstBelow"; threshold = exp(-1)

mdDelay = function(data, nbins = 10, maxlag = 10, criterion = "firstBelow",
                   threshold = exp(-1)){
  
  ## check the data type and possible errors in it
  ## class argument in R 4.0 may return more than a single answer
  ## we just consider the first one as valid
  
  tdata = class(data)[1]
  
  if (tdata == "data.frame"){data = as.matrix(data)} ## convert data.frames into matrices
  
  if (sapply(data, is.numeric)[1] != TRUE){
    stop('Input is not numeric')}
  
  if (is.vector(data) == TRUE & length(data) <= 1){
    stop('Input must be a vector or matrix')}
  
  if (is.matrix(data) == TRUE & ncol(data) <= 1){
    stop('Input must be a vector or matrix')}
  
  
  ncol = ncol(data);
  
  #
  # Calculation of the mutual information as a function of time lag
  #
  
  # Allocate a matrix, where each column will be the auto mutual information
  # as a function of time lag [0; maxlag] for a variable in the input data.
  
  auto_mi = matrix(0, nrow = maxlag + 1, ncol = ncol);
  
  # Allocate a vector to hold the estimated optimal time lag for each
  # dimension.
  lags = rep(0, ncol);
  
  c = 1
  for (c in 1:ncol){
    auto_mi[,c] = autoMI(data[, c], nbins, maxlag);
    if (criterion == 'firstBelow'){
      lags[c] = findFirstBelowThreshold(auto_mi[, c], threshold);  ## double check this function
    } 
    if (criterion == 'localMin'){
      lags[c] = findFirstLocalMinimum(auto_mi[, c]);
    }
  }
  
  delay = mean(lags, na.rm = TRUE);
  
  if (delay == 0){
    stop("Delay not found, consider increasing the threshold, 
         or explore more lags on a larger number of bins")
  }
  
  return(delay = round(delay))
  
}







