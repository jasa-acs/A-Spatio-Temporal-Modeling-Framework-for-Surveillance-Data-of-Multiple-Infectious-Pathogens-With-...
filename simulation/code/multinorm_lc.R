### function for sampling from multivariate distribution with singular precision matrix

library(mvtnorm)

multiNorm <- function(k, mean, P){
n <- nrow(P)
I <-  diag(n-1)
Z <- rmvnorm(k, rep(0, (n-1)), I, method = "chol")
Eigen <- eigen(P)
U <- Eigen$vector[,-n]
Lambda <- diag(1/sqrt(Eigen$values[-n]))
alpha <- U%*%Lambda%*%t(Z)+mean
return(alpha)
}
