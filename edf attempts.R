# Tue Feb 13 17:17:37 2024 ------------------------------


edf <- function(X, lambda, S) {
  edfmat <- MASS::ginv(t(X)%*%X + lambda[1] * S[[1]] + lambda[2] * S[[2]])%*%t(X)%*%X
}  

X <- gamprep$pregam$X[,-1]
lambda <- out$mean$lambda
S <- gamprep$pregam$S


S[[1]] * lambda[1]

myedf <- edf(X, lambda, S)
myedf
sum(myedf)
sum(diag(myedf))
