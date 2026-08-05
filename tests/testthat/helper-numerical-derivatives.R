central_derivative <- function(f, x, h = 1e-6) {
  (f(x + h) - f(x - h)) / (2 * h)
}

central_second_derivative <- function(f, x, h = 1e-4) {
  (f(x + h) - 2 * f(x) + f(x - h)) / h^2
}

central_gradient <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(1))
}

central_hessian <- function(f, x, h = 1e-4) {
  p <- length(x)
  hessian <- matrix(0, p, p)
  f0 <- f(x)
  for (j in seq_len(p)) {
    xp <- x
    xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    hessian[j, j] <- (f(xp) - 2 * f0 + f(xm)) / h^2
    if (j < p) {
      for (k in (j + 1L):p) {
        xpp <- xpm <- xmp <- xmm <- x
        xpp[j] <- x[j] + h
        xpp[k] <- x[k] + h
        xpm[j] <- x[j] + h
        xpm[k] <- x[k] - h
        xmp[j] <- x[j] - h
        xmp[k] <- x[k] + h
        xmm[j] <- x[j] - h
        xmm[k] <- x[k] - h
        hessian[j, k] <- hessian[k, j] <-
          (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h^2)
      }
    }
  }
  hessian
}
