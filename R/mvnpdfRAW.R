#' Multivariate normal distribution density
#'
#' Description
#'
#' @param x a matrix with n columns (the observations) and p rows
#' @param mean a vetor of means
#' @param varcovM a variance-covariance matrix
#' @param Log a logical parameter, with default value to TRUE
#'
#' @import mvtnorm
#' @return a list containing the matrix x, and a vector of length n of the multivariate normal distribution density
#' @export
#' @examples
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}
