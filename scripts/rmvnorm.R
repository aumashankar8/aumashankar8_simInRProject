#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 4: Introduction to Monte Carlo Simulation
#
#' Discuss what this function does...
#'
#' @param n Describe this parameter here...
#' @param mu Describe this parameter here...
#' @param Sigma Describe this parameter here...
#' @return Describe what this returns...
rmvnorm <- function(n, mu, Sigma) {
  if (NROW(mu) != NROW(Sigma)) {
    stop("Matrix Dimensions do not match")
  }
  Z <- matrix(rnorm(n*NROW(mu), 0, 1), ncol = NROW(mu))
  ranData <- matrix(1, n) %*% t(mu) + Z %*% chol(Sigma)
  return(ranData)
}