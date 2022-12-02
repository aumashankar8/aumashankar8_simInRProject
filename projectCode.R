#-------------------------------
#-------------------------------
# Sim in R - Final Project
# Simple HLM and Monte-Carlo
# Akshay Umashankar
# au3692
#-------------------------------
## Libraries ##
library(lme4)
library(rblimp)

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
  
  draw <- rmvnorm(n, c(p_x$mux, 0), covMat) #DRAWS X1 X2 and E
  
  beta0 <- p_y$beta0
  Y_method1 <- beta0 + p_y$betaMat[1] * draw[, 1] + p_y$betaMat[2] * draw[, 2] + draw[, NCOL(draw)]
  
  output <- data.frame(
    X = draw[, -1],
    Y = Y_method1
    #weight = weights
  )
}

montHelp <- monteFunction(n, p_x, p_y)

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

############# Data Generation ################
coefUnconditional <- c()
coefConditional <- c()

parameters <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = c(0.10, 0.25),
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.5
)

dataGenerated <- function(parameters){
  l1_var <- with(parameters, varianceY*(1-ICC))
  l2_varTotal <- with(parameters, varianceY*ICC)
  gamma1 <- with(parameters, d*sqrt(varianceY))
  
  treat <- with(parameters, c(rep(1, J/2), rep(0, J/2))) #Set up 50 Individuals in treatment and 50 individuals No treatment for 1 cluster
  
  L2_id <- sort(with(parameters, rep(c(1:J), nj)))
  L1_id <- with(parameters, rep(c(1:nj), J))
  X <- cbind(L2_id, L1_id, treat)
  Y_within <- with(parameters, rnorm(J*nj, mean, l1_var))
  l2_varExplained <-  0.25*gamma1^2 
  l2_residual <- l2_varTotal - l2_varExplained
  Y_between <- with(parameters, gamma1 * treat + rnorm(J, 0, sqrt(l2_residual)))
  L2_Y <- Y_between %x% rep(1, nj)
  
  Y <- Y_within + L2_Y
  data <- data.frame(X, Y)
  #colnames(data, c("L2_id", "L1_id", "treat", "Y"))
  return(data)
}

dataGenerated(parameters)


############# Data Analysis ################

#Unconditional
uncondMonte <- lmer(Y ~ 1 + (1|L2_id), data)
model1 <- summary(uncondMonte)
coefUnconditional <- c(coefUnconditional, model1$coefficients)


#Conditional Model
condMonte <- lmer(Y ~ 1 + treat + (1|L2_id), data)
model2 <- summary(condMonte)
model2$coefficients
head(data)

############ Replications #################

dataArray <- replicate(dataGenerated(parameters), n = 500)




