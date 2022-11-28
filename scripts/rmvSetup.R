#setup Function if Needed
rmvSetup <- function(muX, sigma, rho){
  sigmaX <- diag(sigma) %*% matrix(c(1, rho, rho, 1),2) %*% diag(sigma)
  setup <- list(
    sigmaX = sigmaX,
    muX = muX
  )
  return(setup)
}