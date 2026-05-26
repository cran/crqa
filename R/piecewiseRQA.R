## piece-wise RQA
# typeRQA = c("full", "diagonal") ## TODO for windowed method

## workers: integer ; number of parallel workers (default: 1 = serial).
##          Set > 1 to enable parallel execution via future/furrr.
##          Each worker computes one (i, j) block independently.

piecewiseRQA <- function(ts1, ts2, blockSize, delay = 1, embed = 1, rescale = 0,
                         radius = 0.001, normalize = 0, mindiagline = 2, minvertline = 2,
                         tw = 0, whiteline = FALSE, recpt = FALSE,
                         side = "both", method = "crqa", metric = "euclidean",
                         datatype = "continuous", typeRQA = "full",
                         windowsize = NA,
                         workers = 1L,
                         rr_denom = "full"){
  rr_denom <- match.arg(rr_denom, c("full", "valid"))

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

  ## initialise all parameters value if missing
  if(exists("embed"))         embed = embed else embed = 1
  if(exists("delay"))         delay = delay else delay = 1
  if(exists("normalize"))     normalize = normalize else normalize = 0
  if(exists("radius"))        radius = radius else radius = 1
  if(exists("mindiagline") & mindiagline > 2) mindiagline = mindiagline else mindiagline <- 2
  if(exists("minvertline") & minvertline > 2) minvertline = minvertline else minvertline <- 2
  if(exists("typeRQA"))     typeRQA = typeRQA else typeRQA = "full"

  ## Apply normalize globally over the full series before block decomposition
  ## so that all blocks share the same scale and the assembled RP is consistent.
  if (normalize > 0) {
    switch(normalize,
      { ts1 <- ts1 - min(ts1); ts1 <- ts1 / max(ts1)
        ts2 <- ts2 - min(ts2); ts2 <- ts2 / max(ts2) },
      { ts1 <- scale(ts1); ts2 <- scale(ts2) }
    )
    normalize <- 0L
  }

  nrows = max(1:length(ts1) - delay*(embed-1))
  ncols = max(1:length(ts2) - delay*(embed-1))

  piecewiseRP = Matrix(0, nrow = nrows, ncol = ncols, sparse = TRUE);

  nrs = seq(blockSize, length(ts1) - delay * (embed-1), blockSize)
  nrc = seq(blockSize, length(ts2) - delay * (embed-1), blockSize)

  ## build the unified list of (i_val, j_val) pairs — main blocks first
  pairs <- vector("list", length(nrs) * length(nrc))
  k <- 1L
  for (ii in nrs) {
    for (jj in nrc) {
      pairs[[k]] <- c(ii, jj)
      k <- k + 1L
    }
  }
  pairs <- pairs[seq_len(k - 1L)]

  ## edge / corner blocks (mirrors original conditional)
  i_last <- if (length(nrs) > 0L) nrs[length(nrs)] else 0L
  j_last <- if (length(nrc) > 0L) nrc[length(nrc)] else 0L
  er = length(ts1) - delay*(embed-1)
  ec = length(ts2) - delay*(embed-1)

  has_edges <- (length(ts1) - i_last - delay*(embed-1) > 0 |
                length(ts2) - j_last - delay*(embed-1) > 0)

  if (has_edges) {
    for (jj in nrc)  pairs <- c(pairs, list(c(er, jj)))
    for (ii in nrs)  pairs <- c(pairs, list(c(ii, ec)))
    pairs <- c(pairs, list(c(er, ec)))
  }

  ## set up parallel workers if requested, restore plan on exit
  if (workers > 1) {
    suppressPackageStartupMessages({
      old_plan <- plan(multisession, workers = workers)
    })
    on.exit(plan(old_plan), add = TRUE)
  }

  ## per-block computation
  .compute_block <- function(pair) {
    i_val <- pair[1]; j_val <- pair[2]
    win1 = ts1[(1 + (i_val - blockSize)):(i_val + delay * (embed - 1))]
    win2 = ts2[(1 + (j_val - blockSize)):(j_val + delay * (embed-1))]
    temp = crqa(win1, win2, delay, embed, rescale,
                radius, normalize, mindiagline, minvertline,
                tw, whiteline, recpt, side, method, metric,
                datatype, rr_denom = rr_denom)$RP
    list(i_val = i_val, j_val = j_val, RP = temp)
  }

  ## execute serially or in parallel
  if (workers > 1) {
    block_results <- future_map(
      pairs, .compute_block,
      .options = furrr_options(seed = TRUE)
    )
  } else {
    block_results <- lapply(pairs, .compute_block)
  }

  ## assemble piecewiseRP (sequential — writes to distinct non-overlapping blocks)
  for (res in block_results) {
    if (!is.null(res$RP)) {
      piecewiseRP[(1 + (res$j_val - blockSize)):res$j_val,
                  (1 + (res$i_val - blockSize)):res$i_val] <- res$RP
    }
  }

  ## extract RQA measures from the assembled piecewise RP
  if (!typeRQA %in% c("full", "diagonal"))
    stop(sprintf('piecewiseRQA(): typeRQA must be "full" or "diagonal", got "%s"', typeRQA))

  if (typeRQA == "full"){
    recpt = TRUE
    results = crqa(piecewiseRP, ts2 = NA, delay, embed, rescale,
                   radius, normalize, mindiagline, minvertline,
                   tw, whiteline, recpt, side, method, metric,
                   datatype, rr_denom = rr_denom)
  }

  if (typeRQA == "diagonal"){
    RP = piecewiseRP

    if (is.logical(RP) != T){
      len      = nrow(RP)-1
      lags     = -len:len

      RP      = matrix(as.numeric(RP), nrow = nrow(RP), ncol = ncol(RP))
      RP_form = unlist(lapply(split(RP, row(RP) - col(RP)), mean, na.rm = TRUE))
      wn      = which(lags >= -windowsize &  lags <= windowsize)

      drpd    = RP_form[wn]
    } else {
      drpd    = 0
    }

    maxrec = max(drpd);
    maxlag = which(drpd == maxrec);

    results = list(profile = drpd, maxrec = maxrec, maxlag = maxlag)
  }

  return(results)

}
