## wincrqa(): Function to run windowed recurrence analyses
##
## Arguments:
##      ts1                     : numeric vector or matrix ; first time series or full recurrence matrix
##      ts2                     : numeric vector ; second time series
##      windowstep              : numeric ; the interval by which the window is moved
##      windowsize              : numeric ; the size of the window
##      delay                   : numeric ; the lag for creating copies in phase-space reconstruction (i.e., lag intervals)
##      embed                   : numeric ; the number of embedding dimensions in phase-space reconstruction (i.e., number of copies)
##      radius                  : numeric ; the threshold distance by which points are considered recurrent;
##                                (default: .001)
##      rescale                 : numeric; choice for rescaling the distance matrix:
##                                    rescale = 0: no rescaling;
##                                    rescale = 1: rescale by mean distance of the entire matrix;
##                                    rescale = 2: rescale by maximum distance of the entire matrix;
##                                    rescale = 3: rescale by minimum distance of the entire matrix;
##                                    rescale = 4: rescale by Euclidean distance of the entire matrix;
##                                (default: 0)
##      normalize               : numeric ; choice for normalizing the time series:
##                                    normalize = 0: no normalizing;
##                                    normalize = 1: normalize by unit interval;
##                                    normalize = 2: normalize by z-score;
##                                (default: 0)
##      mindiagline             : numeric ; minimum number of sequential diagonal points to be considered a line
##                                (default: 2)
##      minvertline             : numeric ; minimum number of sequential vertical points to be considered a line
##                                (default: 2)
##      tw                      : numeric ; the Thieler window parameter, which removes the center diagonal line(s)
##                                (default: 0)
##      whiteline               : boolean ; whether or not to calculate empty vertical lines
##                                (default: FALSE)
##      recpt                   : boolean ; whether the input (ts1) is a recurrence matrix or not
##                                (default: FALSE)
##      side                    : character ; choice of side for the recurrence plot the measures are calculated (relative to center line):
##                                    side = "upper": upper half of the RP (from the center line and excluding the center line)
##                                    side = "lower": lower half of the RP (from the center line and excluding the center line)
##                                    side = "both": entire RP (including the center line)
##                                (default: "both")
##      method                  : character ; choice of recurrence analysis to perform:
##                                    method = "rqa": autorecurrence
##                                    method = "crqa": cross-recurrence
##                                    method = "mdcrqa": multidimensional cross-recurrence
##                                (default: "crqa")
##      metric                  : character ; choice of the type of distance metric used; argument will accept any method from rdist()
##                                (default: "euclidean")
##      datatype                : character ; the kind of data being analyzed:
##                                    datatype = "categorical": categorical data
##                                    datatype = "continuous": continuous data
##                                (default: "continuous")
##      trend                   : boolean ; whether or not to calculate the change in density of recurrence from the center line outward
##                                (default: FALSE)
##      workers                 : integer ; number of parallel workers (default: 1 = serial).
##                                Set > 1 to enable parallel execution via future/furrr.
##                                RAM warning: each worker runs a full crqa() call. Keep workers
##                                low if time series are long (Stage 2 will fix the RAM bottleneck).
## Value:
##     A dataframe of recurrence quantification metrics for each window
##
## Author(s): Moreno I. Coco
## Contributed by: Alexandra Paxton
## Acknowledgements: based on the original MATLAB code by Rick Dale

.packageName <- 'crqa'

wincrqa <- function(ts1,
                    ts2,
                    windowstep,
                    windowsize,
                    delay,
                    embed,
                    radius = 0.001,
                    rescale = 0,
                    normalize = 0,
                    mindiagline = 2,
                    minvertline = 2,
                    tw = 0,
                    whiteline = FALSE,
                    recpt = FALSE,
                    side = 'both',
                    method = 'crqa',
                    metric = 'euclidean',
                    datatype = 'continuous',
                    trend = FALSE,
                    workers = 1L,
                    rr_denom = "full"){
  rr_denom <- match.arg(rr_denom, c("full", "valid"))

  ## stop immediately if the windowsize is smaller than delay*phase AND we're not supplying an RP
  if ((windowsize < embed*delay) & recpt == FALSE){
    stop("Phase-space (embed*delay) longer than windowsize")
  }

  # initialize a flag that will be used to check whether the last window is included
  irregular = FALSE

  ## check the different contexts in which the analyses are running
  if (recpt == FALSE){

    # for RQA: check that the user as provided the same data to allow windowed recurrence
    if (method == "rqa"){
      if (exists("ts2")) {
        ts2 = ts2
      } else {
        stop("Please provide the same vector as argument of t2")
      }
      maxd = length(ts1) ## the total number of points
    }

    # for CRQA
    if (method == "crqa"){
      ts1 = as.vector(as.matrix(ts1))
      ts2 = as.vector(as.matrix(ts2))
      maxd = length(ts1) ## the total number of points

    }

    # for MdCRQA
    if (method == "mdcrqa"){
      ts1 = as.matrix(ts1)
      ts2 = as.matrix(ts2)
      maxd = nrow(ts1) # the total number of points
    }
  } else {
    # for calculating from a supplied RP matrix
    logical_RP = as.logical(ts1)
    v1l = nrow(ts1)
    ts1 = matrix(logical_RP,
                 ncol = as.numeric(v1l))
    v1l = nrow(ts1)
    v2l = ncol(ts1)
    ind = which(ts1 > 0, arr.ind = TRUE)
    r = ind[,1]
    c = ind[,2]

    # total possible number of points on the RP
    maxd = v1l
  }


  ## Apply normalize globally over the full series so every window shares
  ## the same scale. After this block normalize is set to 0 so the per-window
  ## crqa() calls skip the step. (rescale is intentionally left per-window:
  ## it depends on each window's distance matrix, so a global value cannot
  ## be pre-computed without materialising the full N×N matrix.)
  if (normalize > 0 && !recpt) {
    switch(normalize,
      { ## 1: unit interval — mirrors crqa() logic exactly
        ts1 <- ts1 - min(ts1); ts1 <- ts1 / max(ts1)
        ts2 <- ts2 - min(ts2); ts2 <- ts2 / max(ts2) },
      { ## 2: z-score
        ts1 <- scale(ts1)
        ts2 <- scale(ts2) }
    )
    normalize <- 0L
  }

  ## need to make sure to include also "irregular" final segments
  points = seq(1, (maxd - (windowsize)-1), windowstep)

  fpoint = points[length(points)] ## beginning of last window

  ## as the dimension of the data may often be not that precise
  if (fpoint != (fpoint + windowsize)){
    irregular = TRUE ## flag used later on
    fpoint = fpoint + windowsize # update last point
    points = c(points,  fpoint) ## add final point
  }

  ## set up parallel workers if requested, restore plan on exit.
  ## suppressPackageStartupMessages avoids the noisy "package 'future'
  ## was built under R version x.y.z" warning that some R installs emit
  ## the first time future is loaded.
  if (workers > 1) {
    suppressPackageStartupMessages({
      old_plan <- plan(multisession, workers = workers)
    })
    on.exit(plan(old_plan), add = TRUE)
  }

  ## capture internal helpers so the closure is self-contained on workers
  .diags_extract <- diags_extract

  ## per-window computation — runs identically serial or parallel
  .compute_window <- function(idx) {
    i <- points[idx]

    if (irregular == TRUE & i == fpoint) {
      ixs = i:maxd
    } else {
      ixs = i:(i + windowsize - 1)
    }

    if (recpt == FALSE){
      if (method == "crqa" | method == "rqa"){
        ts1win = ts1[ixs]
        ts2win = ts2[ixs]
      }
      if (method == "mdcrqa"){
        ts1win = ts1[ixs, ]
        ts2win = ts2[ixs, ]
      }
    } else {
      ts1win = as.matrix(ts1[ixs, ixs])
      ts2win = NA
    }

    if ((length(ts1win) > delay*embed) || recpt == TRUE){

      ans = crqa(ts1win, ts2win, delay, embed, rescale, radius,
                 normalize, mindiagline, minvertline,
                 tw, whiteline, recpt, side, method, metric,
                 datatype, rr_denom = rr_denom)

      RP = ans$RP
      if (length(RP) == 1) { RP = vector() }

      ans = as.numeric( unlist(ans[1:11]) )

      if (trend == TRUE){

        if (length(RP) > 0){
          RP = as.matrix(RP)
          NX = ncol(RP)
          T = vector("numeric", length = NX-1)

          for (k in 1:(NX-1)){
            ixi_rec = which(.diags_extract(RP, k) != F)
            T[k] = length(ixi_rec)/ (NX-k)*100;
          }

          Ntau = NX - 1 - round(0.1*NX);
          p = pracma::polyfit(2:(Ntau+1), T[1:Ntau], 1)
          TREND = 1000 * p[1]
        } else { TREND = NA }
      } else {
        TREND = NA
      }

      return( c(ans, idx, TREND) )

    } else {
      warning( paste("Window", idx, "was removed because windowsize < delay*embed", sep = " ") )
      return(NULL)
    }
  }

  ## execute serially or in parallel
  if (workers > 1) {
    results_list <- future_map(
      seq_along(points), .compute_window,
      .options = furrr_options(seed = TRUE)
    )
  } else {
    results_list <- lapply(seq_along(points), .compute_window)
  }

  ## drop skipped windows and assemble
  results_list <- Filter(Negate(is.null), results_list)
  crqwin <- do.call(rbind, results_list)

  ## name the measures
  colnames(crqwin) = c("RR", "DET", "NRLINE", "maxL", "L",
                       "ENTR", "rENTR", "LAM", "TT", "catH",
                       "max_vertlength", "win", "TREND")

  return(as.data.frame(crqwin))

}
