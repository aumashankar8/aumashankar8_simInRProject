#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 5: Introduction to Maximum Likelihood
#
#' Discuss what this function does...
#'
#' @param n Describe this parameter here...
#' @param p_x Describe this parameter here...
#' @param p_y Describe this parameter here...
#' @return Describe what this returns...
generate_reg <- function(n, p_x, p_y) {
  #First setup the params for rmvNorm
  rmvSetup <- function(muX, sigma, rho){
    R_mat <- diag(NROW(diag(p_x$sigmax)))
    R_mat[lower.tri(R_mat)] <- rho
    R_mat[upper.tri(R_mat)] <- rho
    sigmaX <- diag(sigma) %*% R_mat %*% diag(sigma)
    setup <- list(
      sigmaX = sigmaX,
      muX = muX
    )
    return(setup)
  }
  setupParam <- with(p_x, rmvSetup(mux, sigmax, rho))
  x_Mat <- rmvnorm(n, setupParam$muX, setupParam$sigmaX)
  
  #Method 1 Generation
  varE <- with(p_y, t(betaMat) %*% setupParam$sigmaX %*% betaMat %*% (1/r_squared - 1))
  beta0 <- with(p_y, muY - t(setupParam$muX) %*% betaMat)
  Y_method1 <- with(p_y, rnorm(n, beta0[1,1] + x_Mat %*% p_y$betaMat, sqrt(varE)))
  output <- ols_regression(Y_method1, cbind(1, x_Mat))
  attr(output, 'parameters') <- c(
    beta0 = beta0,
    betas = p_y$betaMat,
    std = sqrt(varE),
    r_squared = p_y$r_squared
  )
  return(list(
    output = output,
    X = x_Mat,
    Y = Y_method1
  ))

}
