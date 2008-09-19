#############################################################
#
#	quantile.mahalanobis function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 28, 2008
#	Version: 0.2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.mahalanobis <- function(x, probs, all=FALSE) {
  x <- data.matrix(x)
  if (ncol(x) < 2) stop('At least two columns must be supplied')
  nx <- nrow(x)
  S <- cov(x)
  if (abs(det(S)/2) < .Machine$double.eps) stop('The covariance matrix seems to be singular')
  Sinv <- solve(S)
  mah <- rep(0, (nx^2-nx)/2) 
  k <- 0
  for(i in 1:(nx-1)) {
    for(j in (i+1):nx) {
      k <- k+1
      temp <- x[i,]-x[j,]
      mah[k] <- sqrt(temp%*%Sinv%*%temp)
    }
  }
  res <- quantile.default(mah, probs)
  if (all)
    res <- list(quantile=res, stats=mah, call=match.call())
  return(res)
}
