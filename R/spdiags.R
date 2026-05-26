.packageName <- 'crqa'

spdiags <- function(A) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  m <- nrow(A)
  n <- ncol(A)
  k <- min(m, n)

  nz <- which(A != 0, arr.ind = TRUE)

  if (nrow(nz) == 0L) {
    return(list(B = matrix(0.0, nrow = k, ncol = 0L), d = integer(0L)))
  }

  rows     <- nz[, 1L]
  cols     <- nz[, 2L]
  d_elem   <- cols - rows                  # diagonal offset (col - row)
  active_d <- sort(unique(d_elem))
  p        <- length(active_d)

  # Position within each B column follows MATLAB spdiags convention:
  #   fat/square (m <= n): row index of element
  #   tall       (m >  n): col index of element
  pos <- if (m <= n) rows else cols
  grp <- match(d_elem, active_d)           # column in B for each nonzero

  B        <- matrix(0.0, nrow = k, ncol = p)
  B[cbind(pos, grp)] <- A[nz]

  list(B = B, d = active_d)
}
