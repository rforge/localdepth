#############################################################
#
#	localdepth.similarity.mahalanobis function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: October, 02, 2008
#	Version: 0.2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.similarity.mahalanobis <- function(x, y=NULL, tau, weight=NULL) {
  mahdepth <- function(x, mean, covinv) {
    temp <- (x-mean)
    1/(1+temp%*%covinv%*%temp)
  }
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, nrow=1)
  if (is.vector(y))
    y <- matrix(y, nrow=1)
  x <- as.matrix(x)  
  y <- as.matrix(y)
  
  if (ncol(x) < 2) stop('At least two columns must be supplied')
  nx <- nrow(x)
  ny <- nrow(y)
  S <- cov(x)
  media <- as.vector(apply(x, 2, mean))
  if (abs(det(S)/2) < .Machine$double.eps) stop('The covariance matrix seems to be singular')
  Sinv <- solve(S)
  mah <- matrix(0,ny,nx);
  for(i in 1:ny) {
    for(j in 1:nx) { 
      temp <- y[i,]-x[j,]
      mah[i,j] <- sqrt(temp%*%Sinv%*%temp)
    }
  }
  simld <- matrix(0, ny, ny)
  for(i in 1:ny) {
    vetti <- mah[i,]
    for(j in i:ny) {
      vettj <- mah[j,]
      simld[i,j] <- simld[j,i] <- sum(vetti <= tau & vettj <= tau)/nx
    }
  }

  if (!is.null(weight)) {
    mah <- weight(mah)
    simld <- simld*mah
  }
  
  result <- list()
  result$localdepth <- simld
  result$depth <- NA
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- NA
  result$call <- match.call()
  result$tau <- tau
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'mahalanobis'
  class(result) <- 'localdepth.similarity'
  return(result)
}
