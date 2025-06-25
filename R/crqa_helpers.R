## several helpers functions

.packageName <- 'crqa'

###############################################################

ami <- function(ts1, ts2, lag){
  
  lag = round(lag);   ## make sure that lags are integers
  
  ts1 = as.vector(ts1)
  ts2 = as.vector(ts2)
  
  n = length(ts1)
  
  if (n != length(ts2)){
    stop('ami() ts1 and ts2 should have the same length');
  }
  
  ts1 = ts1-min(ts1)
  ts1 = ts1*(1-eps(1))/max(ts1)
  
  ts2 = ts2-min(ts2)
  ts2 = ts2*(1-eps(1))/max(ts2)
  
  v = vector(mode = "numeric", length = length(lag))
  
  lastbins = 0
  
  ii = 2
  
  for (ii in 1:length(lag)){
    
    abslag = abs(lag[ii]);
    
    ## Define the number of bins
    bins = floor(1+log2(n-abslag)+0.5)
    
    if (bins != lastbins){
      bints1 = floor(ts1*bins)+1
      bints2 = floor(ts2*bins)+1
    }
    
    lastbins = bins
    
    PtS = matrix(0, nrow = bins, ncol = bins)
    
    for (jj in 1:(n-abslag)){
      kk = jj + abslag
      
      if (lag[ii] < 0) {
        temp = jj; jj = kk; kk = temp; ## swap
      }
      
      PtS[bints1[kk], bints2[jj]] = PtS[bints1[kk],bints2[jj]] + 1;
    }
    
    PtS = PtS/(n-abslag)
    PtS = PtS + eps(1) ## avoid division and log of zero
    
    Pts1 = rowSums(PtS)
    Pts2 = colSums(PtS)
    
    q = PtS/ outer(Pts1,Pts2);
    q = PtS * log2(q);
    
    v[ii] = sum(q)/log2(bins);
    
  }
  
  return(v)
  
}

# ===============================================================
# autoMI

# used in function mdDelay() to compute the average mutual information
# ===============================================================
# Matlab creator: Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
#     Dan Moenster, Aarhus University
# R translation by: by Moreno I. Coco (moreno.cocoi@gmail.com)

autoMI = function(x, nbins, maxlag){
  
  ami = vector('numeric', length = maxlag+1);
  
  # Compute a histogram of x
  ## reconstruct the sequence of breaks
  seqb = seq(min(x), max(x), length.out = nbins+1);
  
  # extract the counts
  p_i = hist(x, breaks = seqb, plot = F)$counts;
  # p_i = histc(x, edges =  seqb), plot = F)$counts;
  # Normalize to make it a probability
  p_i = p_i / sum(p_i);
  
  # Save the original input vector
  xorig = x;
  
  # Loop over time lags from 0 to maxlag
  tau = 1;
  for (tau in 0:maxlag){
    
    # Make a time delayed version of x
    xtau = circshift(xorig, -tau);
    
    # Remove the last tau elements of this vector and of x
    xtau = xtau[1:(length(xtau)-tau)];
    x    = xorig[1:(length(xorig)-tau)];
    
    # Compute the two-dimensional histogram of x and xtau
    p_ij = hist2d(x, xtau, nbins, show = F)$counts;
    p_ij = p_ij/sum(p_ij)
    
    # Normalize the rows to get the joint probability the x is
    # in bin i and xtau is in bin j.
    # row_sum = sum(p_ij, 2);
    # If any rows sum to zero, replace the sum by 1
    # row_sum(row_sum == 0) = 1;
    # Divide each element with the row sum.
    # p_ij = diag(1 ./ row_sum) * p_ij;
    # Compute the auto mutual information by summing over all
    # combinations of i and j
    i = 1; j = 1
    for (i in 1:nbins){
      for (j in 1:nbins){
        
        if (p_i[i] != 0 & p_i[j] != 0 & p_ij[i,j] != 0){
          ami[tau+1] = ami[tau+1] + p_ij[i,j] * log(p_ij[i,j] / (p_i[i] * p_i[j]))};
      }
    }
  }
  return(ami)
}

# =============================================================================
# catEnt = compute categorical entropy
# Standalone function we can use on the already extracted RP
# (a Sparse Matrix) in order to compute categorical entropy
# using the findBlocks() function.
#
# creator: Giuseppe Leonardi (giuleonardi@gmail.com)
# integrated and modified by Moreno I. Coco (moreno.cocoi@gmail.com)
# =============================================================================
# ind  = the indeces of the recurrence plot where recurrence is observed
# size = the dimension of the recurrence plot

catEnt <- function(ind_c, size){
  
  # print(dim(ind_c))
  points <- findBlocks(ind_c, size)
  
  nbk   <- length(unique(points$x))
  areas <- as.numeric(table(points[points$x%in%1:nbk, 3]))
  
  bkarea <- areas[areas > 1]
  if (length(bkarea) > 0) {
    tabarea <- as.data.frame(table(bkarea))
    parea <- tabarea$Freq / sum(tabarea$Freq)
    catentropy = -sum(parea * log(parea))
  } else { catentropy <- NA}
  
  return(catentropy)
  
}

# Extract diagonals from a matrix.
#
# x A matrix with more than one row AND more than one column.
# which A single numeric that indicates which diagonal to extract. A value of zero extracts the main diagonal, whereas negative values extract diagonals from the upper triangle and positive values extract diagonals from the lower triangle. Diagonals further from the main diagonal have \code{which} values further from zero. If \code{is.null(which)}, then a matrix of diagonal indices for \code{which} is shown.
# incl.labels A single string that indicates whether \code{"row"}, \code{"column"}, or no (\code{"none"}) labels from \code{x} should be returned with the values on the diagonal. Will return numeric values if the labels are all diagonal, otherwise character labels are returned.
# val.name A single string to name the variable that contains the values from the diagonal in the returned data.frame.
# label.name A single string to name the variable that contains the labels in the returned data.frame (see \code{incl.labels})
# return A data.frame with one variable that contains the values from the chosen diagonal of \code{x} and, optionally, a second variable that contains the chosen labels for those values.
#
# Derek H. Ogle, \email{derek@@derekogle.com}, but relied heavily on \url{http://stackoverflow.com/a/27935808/1123933}.
# revised by Moreno I. Coco and changed name to avoid conflict with defunct version in FSA
# Square numeric matrix
# mat1 <- matrix(1:16,nrow=4)
# colnames(mat1) <- LETTERS[1:ncol(mat1)]
# rownames(mat1) <- 1:nrow(mat1)
# mat1
# diags_extract(mat1,which=NULL)
# diags_extract(mat1)
# diags_extract(mat1,which=-1)
# diags_extract(mat1,which=2)
# diags_extract(mat1,incl.labels="row")
# diags_extract(mat1,which=2,incl.labels="row")
# diags_extract(mat1,which=2,incl.labels="col")
# ( tmp <- diags_extract(mat1,which=2,incl.labels="row",val.name="Freq",label.name="age") )
# str(tmp)

diags_extract <- function(x, which = 0,incl.labels = c("none","row","column"),
                          val.name = "value",label.name = "label") {
  ## check if matrix
  if (!is.matrix(x)) stop("'diags_extract' only works with matrices.")
  
  if (nrow(x)==1 | ncol(x)==1) stop("'x' must have more than 1 row and more than 1 column.")
  ## find indices of diagonals for the matrix
  ## idea from http://stackoverflow.com/a/27935808/1123933
  ind <- row(x)-col(x)
  if (is.null(which)) { # nocov start
    ## Simply show the matrix of indices
    cat("Indices matrix corresponding to 'x'.\n")
    rownames(ind) <- rownames(x)
    colnames(ind) <- colnames(x)
    #print(ind)
    #cat("\n") # nocov end
  } else {
    ## extract diagonal from x according to which
    if (which>max(ind) | which<min(ind)) stop("The 'which' diagonal does not exist in 'x'.")
    res <- x[ind==which]
    ## handle adding names
    incl.labels <- match.arg(incl.labels)
    if (incl.labels=="row")
      res2 <- rownames(x)[apply(ind,MARGIN=1,FUN=function(x) any(x==which))]
    else if (incl.labels=="column")
      res2 <- colnames(x)[apply(ind,MARGIN=2,FUN=function(x) any(x==which))]
    else res2 <- NULL
    ## put together as data.frame and return
    if (!is.null(res2)) {
      suppressWarnings(tmp <- as.numeric(res2))
      if (all(!is.na(tmp))) res2 <- tmp
      res <- data.frame(res2,res,stringsAsFactors=FALSE)
      names(res) <- c(label.name,val.name)
    } else {
      res <- data.frame(res)
      names(res) <- c(val.name)
    }
    res
  }
}


# ===============================================================
# findBlocks
# Matrix transformation by @alexis_laz (Stackexchange)
# Assign to non-zero elements of a sparse Matrix
# (i.e. categorical recurrence plot) a numerical code identifying
# its block membership
# ===============================================================
# creator: Giuseppe Leonardi (giuleonardi@gmail.com)
# integrated and modified by Moreno I. Coco (moreno.cocoi@gmail.com) & Alexandra Paxton (alexandra.paxton@uconn.edu)


findBlocks <- function(ind_c, size) {
  
  # get the total number of possible points in the RP
  lt    = nrow(ind_c)
  
  # convert the index matrix to a list
  ind_c = as.list(as.data.frame(ind_c))
  
  # keep track of when
  blocks = list(lastSeenRow = integer(size[1]),
                lastSeenCol = integer(size[2]),
                gr = integer(lt))
  
  # initialize block counter
  ngr = 0 # initialize the counter of blocks
  
  # cycle through each possible point
  for(k in 1:lt) {
    
    # find the next point in the matrix
    kr <- ind_c$row[k]
    kc <- ind_c$col[k]
    #print(c(k, kr, kc))
    i <- blocks$lastSeenRow[kr]
    j <- blocks$lastSeenCol[kc]
    
    # identify which points actually exist and are contiguous
    if (!is.na(i) && i > 0 && (abs(kc - ind_c$col[i]) == 1))      {
      blocks$gr[k] = blocks$gr[i]
    } else if ((!is.na(j)) && j > 0 && ((abs(kr - ind_c$row[j]) == 1))) {
      blocks$gr[k] = blocks$gr[j]
    } else {
      ngr <- ngr + 1L;
      blocks$gr[k] = ngr
    }
    
    blocks$lastSeenRow[kr] <- k
    blocks$lastSeenCol[kc] <- k
  }
  
  return(data.frame(i = ind_c$row, j = ind_c$col, x = blocks$gr))
  
}

# lineprof(findBlocks(ind_c, size))

# ===============================================================
# findFirstBelowThreshold

# used in function mdFnn() to find the first element below the threshold.
# Then test whether an element below the threshold was found,
# and recover if this is not the case.
# ===============================================================
# Matlab creator: Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
#     Dan Moenster, Aarhus University
# R translation by: by Moreno I. Coco (moreno.cocoi@gmail.com)


findFirstBelowThreshold = function(ami, threshold){
  
  idx = which(ami < threshold)[1];
  
  if (length(idx) == 0){
    print('mdDelay() - findFirstBelowThreshold() No value below threshold found.
          Will use local minimum instead');
    
    # If there is more than one elemtent that has the minimum value
    # the min() function returns the first one.
    lag = findFirstLocalMinimum(ami);
    
  } else {
    #  A value of the index idx = 1 corresponds to lag = 0, so 1 is
    #  subtracted from the index to get the lag.
    # MIC: note sure this make sense, at least not in R, as lag = 1 is effectively not lagging anything
    
    lag = idx
    
  }
  
  return(lag)
  
}

# ===============================================================
# findFirstLocalMinimum
# used in function mdFnn to find the local minimum of false-nearest neighbourghors in multi-dimensional time-series
# ===============================================================
# Matlab creator: Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
#     Dan Moenster, Aarhus University
# R translation by: by Moreno I. Coco (moreno.cocoi@gmail.com)

findFirstLocalMinimum = function(ami){
  # Find all local minima
  idx = which(diff(ami) > 0);
  
  if (length(idx) > 0){
    # Select the first local minimum
    idx = idx[1];
  } else {
    print('mdDelay() - findFirstLocalThreshold() No local minimum found. Will use global minimum instead')
    
    print('Consider increasing maxlag or use firstbelow')
    # If there is more than one elemtent that has the minimum value
    # the min() function returns the first one.
    idx = which(ami == min(ami));
  }
  # A value of the index idx = 1 corresponds to lag = 0, so 1 is
  # subtracted from the index to get the lag.
  # MIC: note sure this make sense, at least not in R, as lag = 1 is effectively not lagging anything
  
  lag = idx;
  
}

# ===============================================================
# numerify

## transform categorical series into numerical series
## arguments: two character vectors (or matrices)
## it returns the vectors (or matrices) converted into continuos variable
## creator: Moreno I. Coco (moreno.cocoi@gmail.com)

# ===============================================================
numerify <- function( ts1, ts2 ){
  
  if (is.matrix(ts1) | is.matrix(ts2)){
    nwts1 = nwts2  = matrix(0, nrow = nrow(ts1), ncol = ncol(ts1))
  } else {nwts1 = nwts2 = vector()}
  
  objects = sort(unique(c(ts1, ts2 )))
  ids = seq(1, length(objects), 1)
  
  for ( o in 1:length(objects) ){
    indts1 = which(ts1 == objects[o], arr.ind = T)
    if (length(indts1) > 0) nwts1[indts1] = ids[o]
    indts2 = which(ts2 == objects[o], arr.ind = T)
    if (length(indts2) > 0) nwts2[indts2] = ids[o]
  }
  
  #cbind(ts1, nwts1, ts2, nwts2)
  return( list(nwts1 = nwts1, nwts2 = nwts2  ) )
  
}

# ===============================================================
# theiler

##  used in function crqa() to define a theiler window
## across the diagonal of the recurrence plot to remove recurrent points along the diagonal of the matrix
## written by Moreno I. Coco (moreno.cocoi@gmail.com)
## S = a sparse recurrent matrix
## tw = the theiler window

# S = rbind(c(0,1,1,1,0), c(0,0,1,1,1), c(1,1,1,0,0),c(1,1,0,0,1),c(1,1,0,1,0))

# ===============================================================

theiler <- function(S,tw) {
  
  if (tw > nrow(S)) stop ("crqa(); tw() Theiler window larger
                          than number of diagonals")
  
  if (tw > 0){ ## remove the theilers' points
    theilers = vector()
    
    for (t in 1:tw){
      
      if (t == 1){ ## this is the diagonal
        theilers = rbind(theilers, which(row(S) == col(S), arr.ind = TRUE),
                         deparse.level = 0)
      } else {
        
        ## above diagonal
        theilers = rbind(theilers, which(row(S) == col(S) + (t -1), arr.ind = TRUE),
                         deparse.level = 0)
        
        ## below diagonal
        theilers = rbind(theilers, which(row(S) == col(S) - (t -1), arr.ind = TRUE),
                         deparse.level = 0)
        
      }
    }
    
    S[theilers] = 0
    return (S)
  } else { ## leave the matrix untouched
    return(S)
  }
  
}

## ===============================================================
## tt

## extract vertical lines from a recurrence plot,
## calculate laminarity and trapping time.
## part of this code was borrowed from the original
## tt.m function of the crptool written by Norbert Marwan

## build a random x matrix
# r = 100; c = 100
# x = round(matrix(runif(r*c), r, c))
# whiteline = FALSE
# minvertline = 2
# ans = tt(x, minvertline, whiteline)
# x = rbind(c(1,0,1,1,1,1), c(1,1,0,0,1,1), c(0,1,0,1,0,1))

## ===============================================================

tt <- function(x, minvertline, whiteline){
  # require(Matrix)
  
  ##########################################
  ## for black vertical lines
  xb = x ## just copy over x in .m use of double
  ## pad the list row with zeros
  
  xb = rbind(xb, rep(0, ncol(x)), deparse.level = 0) #xb(end+1,:) = 0;
  xb = as.vector(xb) ## make it a vector
  z = diff(xb)
  z0 = which(z == 1)  ## begin of black sequence
  z1 = which(z == -1) ## end of black sequence
  
  ## measure the length of black lines
  ## do some padding on the series
  
  if (z0[1] > z1[1]){
    z0 = c(0, z0)  # add one zero on the top = z0(1:end); z0(1,1) = 0
  }
  
  if (length(z0) > length(z1)){
    z0 = z0[-length(z0)] #(end)=[];
  }
  
  t = sort(z1-z0);
  t1 = t[t >= minvertline]
  
  TT = mean(t1) ## trapping time
  max_vertlength = max(t1)
  
  ## calculate laminarity: the amount of laminar phases in the system.
  nrbln = sum(x) ## nr. of recurrent points
  
  ## total number of vertical lines
  ## including those below the vertline threshold.
  mnvert = which(t1 < minvertline)
  
  if (length(mnvert) != 0){ ## there are vertical lines smaller than minimum
    t1tr = t1[-mnvert]
  } else {
    t1tr = t1
  }
  
  ## else just return t1tr
  
  if(length(t1tr) > 0){ ## there are vertical lines
    lam = (sum(t1tr)/nrbln)*100
  } else {
    lam = 0 }
  
  
  ####################################
  ## for white vertical lines
  
  if (whiteline == TRUE){
    
    xw = as.matrix(x) ## copy over x another time
    ind = which(xw > 0, arr.ind = TRUE) ## extract the indeces with rec. points
    
    r = ind[,1]; c = ind[,2]
    
    vertline = tapply(r, c,
                      function(x){
                        i1 = min(x); i2 = max(x)
                        ind = cbind(i1, i2)
                        return( ind ) } ) ## take min and max indeces of the matrix
    ## where recurrent points start(min) and end(max)
    
    ## as rec. indeces are repeated across columns address them using
    ## vector operations
    
    colind = as.numeric( names(vertline) )
    matind = matrix(as.numeric( unlist(vertline) ),
                    byrow = TRUE, ncol = 2)
    ## put the indeces in a matrix
    matind = rbind(matind, rep(0,ncol(matind))) ## add one column for padding
    reprw = which(diff(matind)[,1] != 0) ## get the length of column
    unind = matind[reprw, ]
    
    
    for (rp in 1:length(reprw) ){
      
      if (rp == 1){
        
        xw[1:unind[rp,1], colind[1:reprw[rp]] ] = 1
        xw[unind[rp,2]:nrow(xw), colind[1:reprw[rp]] ] = 1
        
      } else {
        
        xw[1:unind[rp,1], colind[reprw[rp -1]:reprw[rp] ] ] =  1
        xw[unind[rp,2]:nrow(xw), colind[reprw[rp -1]:reprw[rp] ] ] = 1
        
      }
    }
    
    xw = rbind(xw, rep(1, ncol(xw)), deparse.level = 0)
    
    zw = diff(as.vector(xw)) ;
    
    z0w = which(zw == -1); # begin of white sequence
    z1w = which(zw == 1);  # end of white sequence
    
    
    ## measure the length of white lines
    if (z0w[1] > z1w[1]){
      z0w = z0w[-1]
      if (length(z1w) > length(z0w) ){
        z1w = z1w[-length(z1w)]
      }
    }
    
    
    if ( length(z1w) > length(z0w) ){
      z0w = c(1,z0w)
    }
    
    tw = sort(z1w-z0w)
    t1w = tw[which(tw-1 > 0)]
    
    tw = tw
    
  } else {  tw = NA }
  
  return(list (TT = TT, lam = lam, max_vertlength = max_vertlength, tw = tw, tb = t) )
  
}
