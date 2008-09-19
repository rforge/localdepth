#############################################################
#
#	plot.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 2, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

plot.localdepth <- function(x, xlab="Depth", ylab="Local Depth", main="DD plot", ...) {

  sdepth <- x$depth
  ldepth <- x$localdepth
  sdepth <- (sdepth-min(sdepth))/(max(sdepth)-min(sdepth))
  ldepth <- (ldepth-min(ldepth))/(max(ldepth)-min(ldepth))

  plot(x=sdepth, y=ldepth, xlab=xlab, ylab=ylab, main=main, ...)
}
