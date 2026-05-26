## Approximative Recurrence Quantification Analysis
## See man/aRQA.Rd for full documentation.
##
## Schultz, Spiegel, Marwan & Albayrak (2015) Phys. Lett. A 379:997-1011.
## Spiegel, Schultz & Marwan (2016) in Webber/Ioana/Marwan (eds), Springer.

.packageName <- 'crqa'

aRQA <- function(ts1, ts2 = ts1,
                 delay = 1, embed = 1, radius = 0.1,
                 mindiagline = 2, normalize = 0) {

  if (mindiagline < 1) mindiagline <- 2L
  mindiagline <- as.integer(mindiagline)

  to_mat <- function(x) {
    if (is.matrix(x)) x else matrix(as.numeric(x), ncol = 1L)
  }
  m1 <- to_mat(ts1); m2 <- to_mat(ts2)
  if (ncol(m1) != ncol(m2))
    stop("aRQA(): ts1 and ts2 must have the same number of columns")

  if (nrow(m1) < 1000L)
    warning("aRQA() called on a short series (N = ", nrow(m1), " < 1000). ",
            "method = 'rqa' or 'crqa' will give exact results at comparable ",
            "speed. aRQA is designed for N >> 10,000.", call. = FALSE)
  d <- ncol(m1)

  if (normalize == 1L) {
    m1 <- apply(m1, 2, function(v) { v <- v - min(v); v / max(v) })
    m2 <- apply(m2, 2, function(v) { v <- v - min(v); v / max(v) })
  } else if (normalize == 2L) {
    m1 <- scale(m1); m2 <- scale(m2)
  }

  embed_ts <- function(M, embed, delay) {
    N <- nrow(M) - (embed - 1L) * delay
    if (N <= 0L) stop("aRQA(): time series shorter than (embed-1)*delay")
    if (embed == 1L) return(M[seq_len(N), , drop = FALSE])
    out <- matrix(0, nrow = N, ncol = d * embed)
    for (k in seq_len(embed)) {
      cols <- ((k - 1L) * d + 1L):(k * d)
      out[, cols] <- M[((k-1L)*delay+1L):((k-1L)*delay+N), , drop = FALSE]
    }
    out
  }
  e1 <- embed_ts(m1, embed, delay)
  e2 <- embed_ts(m2, embed, delay)

  RR_L <- function(L) {
    n1 <- nrow(e1) - L + 1L; n2 <- nrow(e2) - L + 1L
    if (n1 <= 0L || n2 <= 0L) return(0)
    if (L == 1L) { s1 <- e1; s2 <- e2 } else {
      cols_each <- ncol(e1); D <- cols_each * L
      s1 <- matrix(0, n1, D); s2 <- matrix(0, n2, D)
      for (l in seq_len(L)) {
        cols <- ((l-1L)*cols_each+1L):(l*cols_each)
        s1[, cols] <- e1[l:(n1+l-1L), , drop = FALSE]
        s2[, cols] <- e2[l:(n2+l-1L), , drop = FALSE]
      }
    }
    ix1 <- floor(s1/radius); ix2 <- floor(s2/radius)
    mins  <- pmin(apply(ix1,2,min), apply(ix2,2,min))
    ix1   <- t(t(ix1)-mins); ix2 <- t(t(ix2)-mins)
    bases <- pmax(apply(ix1,2,max), apply(ix2,2,max)) + 1L
    if (sum(log2(pmax(bases,1))) > 30 || ncol(ix1) > 12) {
      key1 <- do.call(paste, c(as.data.frame(ix1), sep="_"))
      key2 <- do.call(paste, c(as.data.frame(ix2), sep="_"))
      h1 <- table(key1); h2 <- table(key2)
      common <- intersect(names(h1), names(h2))
      same_bin <- sum(as.numeric(h1[common]) * as.numeric(h2[common]))
    } else {
      key1 <- ix1[,1L]; key2 <- ix2[,1L]
      if (ncol(ix1) > 1L) {
        for (k in 2:ncol(ix1)) { key1 <- key1*bases[k]+ix1[,k]; key2 <- key2*bases[k]+ix2[,k] }
      }
      maxkey <- as.integer(prod(bases))
      h1 <- tabulate(key1+1L, nbins=maxkey); h2 <- tabulate(key2+1L, nbins=maxkey)
      same_bin <- sum(as.numeric(h1)*as.numeric(h2))
    }
    same_bin / (as.numeric(n1)*as.numeric(n2))
  }

  RR1    <- RR_L(1L)
  RR_Lm  <- RR_L(mindiagline)
  RR_Lm1 <- RR_L(mindiagline + 1L)

  DET    <- if (RR1 > 0) (mindiagline*RR_Lm - (mindiagline-1L)*RR_Lm1) / RR1 else NA_real_
  dL     <- RR_Lm - RR_Lm1
  Lmean  <- if (!is.na(dL) && dL > 0) (mindiagline*RR_Lm - (mindiagline-1L)*RR_Lm1)/dL else 0

  list(RR=RR1*100, DET=DET*100, NRLINE=NA_integer_, maxL=NA_real_, L=Lmean,
       ENTR=NA_real_, rENTR=NA_real_, LAM=NA_real_, TT=NA_real_,
       catH=NA_real_, max_vertlength=NA_real_, RP=NA)
}
