#############################################################
#
#	quantile.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 07, 2008
#	Version: 0.1-2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.localdepth <- function(x, probs, use=c('volume', 'diameter'),  method=c('simplicial', 'ellipsoid', 'mahalanobis'), nsamp='all', size=FALSE, dimension=NULL, ...) {
  use <- match.arg(use)
  method <- match.arg(method)
## simplicial
  if (method=='simplicial') {
    if (is.circular(x))
      quantile.simp.circular(x=x, probs=probs, all=size)
    else
      quantile.simp(x=x, probs=probs, use=use, nsamp=nsamp, all=size)
## mahalanobis    
  } else if (method=='mahalanobis') {
      if (is.circular(x)) stop("method 'mahalanobis' is not implemented for circular data")
      quantile.mahalanobis(x=x, probs=probs, all=size)
## ellipsoid
  } else {
      if (is.circular(x)) stop("method 'ellipsoid' is not implemented for circular data")
      quantile.ellipsoid(x=x, probs=probs, use=use, nsamp=nsamp, all=size, dimension=dimension)
  }
}
