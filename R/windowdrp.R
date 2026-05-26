## GNU License, written by Moreno I. Coco
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: windowstep = interval considered on the serie;
##            window_size = the size of the window wherin crqa is runned.
##            lag_profile_width = lags within the window

## windowstep = 10;
## windowsize = 100;
## lagwidth = 20;
## datatype = "categorical"

## workers: integer ; number of parallel workers (default: 1 = serial).
##          Set > 1 to enable parallel execution via future/furrr.

.packageName <- 'crqa'

windowdrp <- function(ts1, ts2, windowstep, windowsize, lagwidth,
                      radius = 0.001, delay = 1, embed = 1, rescale = 0,
                      normalize = 0, mindiagline = 2, minvertline = 2,
                      tw = 0, whiteline = FALSE, recpt = FALSE, side = 'both',
                      method = 'crqa', metric = 'euclidean',
                      datatype = 'continuous',
                      workers = 1L,
                      rr_denom = "full"){
  rr_denom <- match.arg(rr_denom, c("full", "valid"))

  irregular = FALSE; # initialize a flag that will be used to check whether the last window is included
  ## check the different contexts in which the analyses are running

  if (method == "rqa"){
    ## check that the user as provided the same data to allow windowed recurrence
    if (exists("ts2")) ts2 = ts2 else stop("Please provide the same vector as argument of t2")
  }

  if (method == "crqa"){
    ts1 = as.vector(as.matrix(ts1));   ts2 = as.vector(as.matrix(ts2))
    maxd = length(ts1) ## the total number of points
    points = seq(1, (maxd - (windowsize)-1), windowstep)
  }

  if (method == "mdcrqa"){
    ts1 = as.matrix(ts1);   ts2 = as.matrix(ts2)
    maxd = nrow(ts1) ## the total number of points
    points = seq(1, (maxd - (windowsize)-1), windowstep)
  }

  ## Apply normalize globally over the full series (see wincrqa for rationale).
  if (normalize > 0) {
    switch(normalize,
      { ts1 <- ts1 - min(ts1); ts1 <- ts1 / max(ts1)
        ts2 <- ts2 - min(ts2); ts2 <- ts2 / max(ts2) },
      { ts1 <- scale(ts1); ts2 <- scale(ts2) }
    )
    normalize <- 0L
  }

  ## need to make sure to include also "irregular" final segments
  fpoint = points[length(points)] ## beginning of last window
  ## as the dimension of the data often may not be that precise
  if (fpoint != (fpoint + windowsize)){
    irregular = TRUE ## flag used later on
    fpoint = fpoint + windowsize # update last point
    points = c(points,  fpoint) ## add final point
  }

  ## set up parallel workers if requested, restore plan on exit
  if (workers > 1) {
    suppressPackageStartupMessages({
      old_plan <- plan(multisession, workers = workers)
    })
    on.exit(plan(old_plan), add = TRUE)
  }

  ## per-window computation
  .compute_window <- function(i) {
    if (irregular == TRUE & i == fpoint) {
      ixs = i:maxd
    } else {
      ixs = i:(i + windowsize - 1)
    }

    if (method == "crqa" | method == "rqa"){
      ts1win = ts1[ixs]
      ts2win = ts2[ixs]
    }

    if (method == "mdcrqa"){
      ts1win = ts1[ixs, ]
      ts2win = ts2[ixs, ]
    }

    mean(drpfromts(ts1win, ts2win, lagwidth,
                   radius, delay, embed, rescale,
                   normalize, mindiagline, minvertline,
                   tw, whiteline, recpt, side,
                   method, metric,
                   datatype, rr_denom = rr_denom)$profile)
  }

  ## execute serially or in parallel
  if (workers > 1) {
    drpd <- unlist(future_map(
      points, .compute_window,
      .options = furrr_options(seed = TRUE)
    ))
  } else {
    drpd <- vapply(points, .compute_window, numeric(1))
  }

  maxrec = max(drpd);
  maxlag = which(drpd == maxrec);

  return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )

}
