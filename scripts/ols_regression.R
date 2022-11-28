#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 5: Introduction to Maximum Likelihood
#
#' Discuss what this function does...
#'
#' @param y Describe this parameter here...
#' @param X Describe this parameter here...
#' @return Describe what this returns...
ols_regression <- function(y, X) {
  beta <- solve(t(X) %*% X) %*% t(X) %*% y #Beta Matrix
  y_hat <- X %*% beta
  #H <- X %*% (solve(t(X) %*% X)) %*% t(X) #hat Matrix
  e <- y - y_hat #residuals
  sdE <- sqrt((t(e) %*% e)/(NROW(y) - (NCOL(X)-1) - 1))
  SE <- sqrt(diag((sdE[1]^2 * solve(t(X) %*% X)))) #Standard Errors for each estimate (diagonal of this matrix)
  R_squared <- (t(beta) %*% cov(X) %*% beta) / (t(beta) %*% cov(X) %*% beta + sdE^2)
  t_val <- beta / SE
  pr <- 2 * pt(abs(t_val), NROW(y) - NCOL(X) - 1, lower.tail = FALSE)
  estimate <- rbind(beta, sdE, R_squared)
  
  #Setup lengths of vectors
  lengths <- max(c(length(estimate), length(SE), length(t_val), length(pr)))
  length(estimate) <- lengths
  length(SE)       <- lengths
  length(t_val)    <- lengths
  length(pr)       <- lengths
  
  
  mat <- cbind(estimate, SE, t_val, pr)
  rownames(mat) <- c(paste0("b", 0:(NROW(mat)-3)), "SD(e)", "R2")
  colnames(mat) <-(c("Estimate", "SE", "t_val","pr"))
  return(mat)
}