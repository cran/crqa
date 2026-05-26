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
                 method = "rqa", metric = "euclidean", datatype = "continuous",
                 rr_denom = "full"){
  
  rr_denom <- match.arg(rr_denom, c("full", "valid"))

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
    chkmet = method%in%c("rqa", "crqa", "mdcrqa", "aRQA")
    if (chkmet == F) stop("The method you have used is not valid")

    ## --- aRQA short-circuit ---------------------------------------------
    ## Stage 3a: approximative RQA via phase-space coarse-graining
    ## (Schultz et al. 2015; Spiegel et al. 2016). Skips the full O(N^2)
    ## pipeline entirely. Parameters tw, side, whiteline, minvertline,
    ## metric, recpt are silently ignored on this path: aRQA approximates
    ## DET/L from histogram statistics rather than scanning a recurrence
    ## matrix. Vertical-line measures (LAM, TT) and full line-length
    ## counts (NRLINE, maxL, ENTR) are returned as NA pending a future
    ## extension; users wanting them should fall back to method="rqa"/"crqa".
    if (method == "aRQA") {
      return(aRQA(ts1 = ts1, ts2 = ts2,
                  delay = delay, embed = embed,
                  radius = radius, mindiagline = mindiagline,
                  normalize = normalize))
    }
    ## ---------------------------------------------------------------------
    
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
    
    ## ---- Fast path: fused Rcpp inner loop -------------------------------
    ## When the metric is one we have a tight C++ kernel for, we skip the
    ## entire cdist -> dm -> threshold -> sparseMatrix -> t() -> theiler ->
    ## side-mask -> line_stats chain and compute distance + threshold +
    ## line statistics in a single O(N*M) pass. Memory drops from O(N^2)
    ## (peak ~2.5x of 8N^2 measured empirically on v2.1.0) to O(N + nnz),
    ## and runtime drops 3-10x at typical sizes.
    fused_metric <- match(metric, c("euclidean", "maximum", "manhattan"))
    use_fused    <- !is.na(fused_metric)

    if (use_fused) {
      ts1_mat <- if (is.matrix(ts1)) ts1 else matrix(as.double(ts1), ncol = 1L)
      ts2_mat <- if (is.matrix(ts2)) ts2 else matrix(as.double(ts2), ncol = 1L)
      storage.mode(ts1_mat) <- "double"
      storage.mode(ts2_mat) <- "double"

      if (rescale > 0) {
        rescaledist     <- crqa_rescale_stat(ts1_mat, ts2_mat,
                                             as.integer(fused_metric),
                                             as.integer(rescale))
        effective_radius <- radius * rescaledist
      } else {
        effective_radius <- radius
      }
      side_id <- match(side, c("both", "upper", "lower")) - 1L
      if (is.na(side_id)) stop("Unknown side: ", side)

      ls_out  <- crqa_fused_cpp(ts1_mat, ts2_mat,
                                 as.double(effective_radius),
                                 as.integer(fused_metric),
                                 as.integer(tw), as.integer(side_id),
                                 as.integer(mindiagline),
                                 as.integer(minvertline),
                                 isTRUE(whiteline))
      ## ls_out$v2l = n1 (pre-transpose nrow), ls_out$v1l = n2 (pre-transpose ncol)
      v1l <- ls_out$v2l
      v2l <- ls_out$v1l
      ## Pre-transpose recurrent indices for catEnt (which expects pre-transpose).
      ## Column names "row"/"col" matter: catEnt/findBlocks index ind by name.
      ind <- if (ls_out$numrecurs > 0L) {
        cbind(row = ls_out$jj, col = ls_out$ii)
      } else {
        m <- matrix(integer(0), 0L, 2L); colnames(m) <- c("row", "col"); m
      }
      r   <- ind[, 1]; c <- ind[, 2]
      used_fused <- TRUE
    } else {
      ## Fallback for metrics not implemented in C++ (canberra, minkowski, ...)
      dm <- as.matrix(cdist(ts1, ts2, metric = metric))
      if (rescale > 0) {
        rescaledist <- switch(rescale,
          mean(dm),
          max(dm),
          min(dm),
          abs(sum(dm) / (nrow(dm) ^ 2 - nrow(dm)))
        )
        effective_radius <- radius * rescaledist
      } else {
        effective_radius <- radius
      }
      v1l <- nrow(dm); v2l <- ncol(dm)
      ind <- which(dm <= effective_radius, arr.ind = TRUE)
      rm(dm)
      r <- ind[, 1]; c <- ind[, 2]
      used_fused <- FALSE
    }
    
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
    used_fused <- FALSE
  }

  if (length(r) != 0 & length(c) != 0){ ##avoid cases with no recurrence

    if (isTRUE(used_fused)) {
      ## Fused path already computed line statistics. Build S from the
      ## recurrent indices for backwards compatibility: $RP has always been
      ## returned regardless of `recpt` (which controls input type, not
      ## output). S construction here is O(nnz log nnz), not O(N^2), so the
      ## RAM win is preserved relative to the legacy path.
      S <- sparseMatrix(i = ls_out$ii, j = ls_out$jj,
                        dims = c(ls_out$v1l, ls_out$v2l))
      numrecurs <- ls_out$numrecurs
      diaglines <- ls_out$diaglines
    } else {
      ## Legacy path (RP-input branch, or fallback metric).
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

      ####################################################################
      ## Stage 2 (2026-05): single sparse-index scan in line_stats().
      ####################################################################
      ls_out = line_stats(S, mindiagline, minvertline, whiteline = whiteline)
      numrecurs = ls_out$numrecurs
      diaglines = ls_out$diaglines
    }

    ## Recurrence-rate denominator selection (rr_denom):
    ##
    ## The literature does not settle on a single convention for the RR
    ## denominator when a Theiler window or one-sided mask is applied.
    ## Both choices below appear in published RQA software and analyses:
    ##
    ##   "full"  (= behaviour of crqa <= 2.0.7, of the original Marwan
    ##            (2007) Phys. Rep. 438 definition  RR = (1/N^2) Sum R_ij,
    ##            and of the TOCSY MATLAB toolbox): the denominator is
    ##            v1l * v2l (all cells), regardless of how many were
    ##            blanked by the Theiler window or the side mask. Those
    ##            blanked cells contribute 0 to the numerator, so RR is
    ##            "diluted" by the exclusion. Useful for reproducing
    ##            published numbers and for apples-to-apples comparison
    ##            with other tools that use this convention.
    ##
    ##   "valid" (= new default in v2.1.0; introduced by pjbruna's
    ##            community PR; generalised here to rectangular RPs via
    ##            theiler_exclusion(v1l, v2l, tw)): the denominator counts
    ##            only the cells that survived the Theiler band and the
    ##            side mask. RR thus reflects the fraction of *analysable*
    ##            cell pairs that are recurrent. Internally consistent
    ##            with how DET/ENTR are already computed in crqa() (their
    ##            denominator is numrecurs, the post-Theiler count of
    ##            recurrent points, not v1l*v2l).
    ##
    ## For tw == 0 and side == "both", both conventions give the same
    ## value. The choice only matters when tw > 0 or side != "both". We
    ## default to "valid" because it is internally consistent with the
    ## DET/ENTR denominators, but the option is provided to recover
    ## historical numbers exactly: pass rr_denom = "full".

    ## Cast to double before multiplication: v1l * v2l overflows int when
    ## v1l*v2l > .Machine$integer.max ~ 2.1e9, i.e. around N >= 46341 for
    ## square RPs.  Without this, RR becomes NA at large N (silent until
    ## downstream code does an `if (region > 0)` test and errors with
    ## "missing value where TRUE/FALSE needed").
    v1l_d <- as.double(v1l)
    v2l_d <- as.double(v2l)
    if (rr_denom == "full") {
      region = v1l_d * v2l_d
    } else if (side == "both") {
      region = v1l_d * v2l_d - as.double(theiler_exclusion(v1l, v2l, tw))
    } else {
      ## "upper": keep cells with j - i >= max(tw, 1);  "lower": i - j >= max(tw, 1)
      min_d = max(as.integer(tw), 1L)
      if (side == "upper") {
        if (min_d > v2l - 1L) {
          region = 0
        } else {
          d_vals = seq.int(min_d, v2l - 1L)
          region = sum(as.double(pmin(v1l, v2l - d_vals)))
        }
      } else { ## "lower"
        if (min_d > v1l - 1L) {
          region = 0
        } else {
          d_vals = seq.int(min_d, v1l - 1L)
          region = sum(as.double(pmin(v2l, v1l - d_vals)))
        }
      }
    }

    percentrecurs = if (region > 0) (numrecurs / region) * 100 else NA_real_

    ## exclude the rare cases where there are no diaglines
    if (length(diaglines) != 0) {

      numdiaglines = length(diaglines)
      maxline = max(diaglines)
      meanline = mean(diaglines)

      tabled = as.data.frame(table(diaglines))

      total = sum(tabled$Freq)
      p = tabled$Freq/total

      ## remove zero probability..it should not be necessary
      del = which(p == 0)
      if (length(del) > 0) {
        p = p[-del]
      }

      ## entropy log2, and relative entropy divided by max
      entropy = - sum(p*log(p))
      relEntropy = entropy/(-1*log(1/nrow(tabled)))

      ## entropy/max entropy: comparable across contexts and conditions.

      pdeter = sum(diaglines)/numrecurs*100
      ## percent determinism: the predictability of the dynamical system

      ## laminarity, trapping time, max vertical line — from line_stats()
      lam = ls_out$lam
      TT  = ls_out$TT
      max_vertlength = ls_out$max_vertlength

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
    
    ## white-line outputs are NA unless whiteline=TRUE was requested.
    ## Three measures are returned:
    ##   wmean = mean white vertical line length
    ##           = approximate mean inter-recurrence time
    ##   wmax  = longest white vertical line
    ##   wENTR = Shannon entropy of the white-line length distribution
    ##           = recurrence-time entropy (RTE), Faure & Korn 2004
    ##           Phys. Lett. A 332:329-339; Little, McSharry, Roberts,
    ##           Costello & Moroz 2007 Biomed. Eng. Online 6:23.
    if (exists("ls_out") && isTRUE(whiteline)) {
      wmean <- ls_out$wmean; wmax <- ls_out$wmax; wENTR <- ls_out$wENTR
    } else {
      wmean <- NA_real_; wmax <- NA_real_; wENTR <- NA_real_
    }

    results = list(RR = percentrecurs, DET = pdeter,
                   NRLINE = numdiaglines, maxL = maxline,
                   L = meanline, ENTR = entropy,
                   rENTR = relEntropy,
                   LAM = lam, TT = TT, catH = catH,
                   max_vertlength = max_vertlength,
                   wmean = wmean, wmax = wmax, wENTR = wENTR,
                   RP = S)

  } else { # print (paste ("No recurrence found") )
    results = list(RR = 0, DET = NA, NRLINE = 0,
                   maxL = 0, L = 0, ENTR = NA, rENTR = NA,
                   LAM = NA, TT = NA, catH = NA,
                   max_vertlength = NA,
                   wmean = NA, wmax = NA, wENTR = NA,
                   RP = NA)}  
  
  return (results)
  
}




