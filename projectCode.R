#-------------------------------
#-------------------------------
# Sim in R - Final Project
# Simple HLM and Monte-Carlo
# Akshay Umashankar
# au3692
#-------------------------------
#-------------------------------
source("scripts/rmvnorm.R")
source("scripts/ols_regression.R")
source("scripts/generate_reg.R")


#Lab 3 Reference
set.seed(88888)
n = 100
p_x <- list(
  rho = 0.3,
  mux = c(0, 0),
  sigmax = c(1, 1),
  rhoxE = 0.3
)

p_y <- list(
  betaMat = c(1, 1),
  beta0 = 100,
  muy = 0,
  r_squared = 0.4
)

monteFunction <- function(n, p_x, p_y){
  rmvSetupMonte <- function(muX, sigma, rho){
    R_mat <- diag(NROW(diag(p_x$sigmax)))
    R_mat[lower.tri(R_mat)] <- rho
    R_mat[upper.tri(R_mat)] <- rho
    sigmaX <- diag(p_x$sigmax) %*% R_mat %*% diag(p_x$sigmax)
    setup <- list(
      sigmaX = sigmaX,
      muX = muX
    )
    return(setup)
  }
  
  setupParam <- with(p_x, rmvSetupMonte(mux, sigmax, rho))
  varE <- with(p_y, t(betaMat) %*% setupParam$sigmaX %*% betaMat %*% (1/r_squared - 1))
  sigmaE <- sqrt(varE)
  
  corMat <- rbind(c(p_x$rhoxE, rep(0, NROW(setupParam$sigmaX) - 1)), setupParam$sigmaX)
  corMat <- cbind(c(1, p_x$rhoxE, rep(0, NCOL(setupParam$sigmaX) - 1)), corMat)
  
  sigmaEMat <- c(sigmaE, rep(1, NCOL(corMat) - 1))
  
  covMat <- diag(sigmaEMat) %*% corMat %*% diag(sigmaEMat)
  
  draw <- rmvnorm(n, c(p_2bxMonte$mux, 0), covMat) #DRAWS X1 X2 and E
  
  beta0 <- p_y$beta0
  Y_method1 <- beta0 + p_y$betaMat[1] * draw[, 1] + p_y$betaMat[2] * draw[, 2] + draw[, NCOL(draw)]
  
  output <- data.frame(
    X = draw[, -1],
    Y = Y_method1
    #weight = weights
  )
}

montHelp <- monteFunction(n, p_2bxMonte, p_2byMonte)

#------------------------------------------
#MSE and Jackknife Functions
factoryMSE <- function(pop_value){
  force(pop_value)
  function(parameter){
    rowMeans((parameter - pop_value)^2)
  }
}

jackknife <- function(data, func){
  i <- 1
  biasVec <- c()
  for (x in 1:length(data[1, 1, ])){
    dataSub <- data[-i]
    biasVec <- c(biasVec, mean(sapply(dataSub, func)))
    i <- i + 1
  }
  thetaJack <- mean(biasVec)
  varJack <- ((length(data)-1)/length(data)) * sum((biasVec - thetaJack)^2)
  sdJack <- sqrt(varJack)
  return(sdJack)
}
#------------------------------------------


#Data Generation
#Level 1 (observations)

#Randomly generate X




#Level 2 (groups)

