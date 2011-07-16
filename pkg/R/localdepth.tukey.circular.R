#############################################################
#
#	localdepth.tukey.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: July, 16, 2011
#	Version: 0.1
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.tukey.circular <- function(x, y=NULL, tau) {
  warning('For now, only the depth is implemented')
  if (is.null(y))
    y <- x
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(y, "class") <- attr(y, "circularp") <-  NULL
  if (!is.vector(x)) stop('y must be a vector')
  if (!is.vector(y)) stop('y must be a vector')
  nx <- length(x)
  if (nx < 2) stop('x must have at least length 2')
  ny <- length(y)
  
  ## number of couples
  
  res <- .Fortran("ldtc",
    x = as.double(x),
    y = as.double(y),
    nx = as.integer(nx),
    ny = as.integer(ny),
    tau = as.double(tau),
    depth = double(ny),
    localdepth = double(ny),
    PACKAGE = "localdepth")
  
  result <- list()
  result$localdepth <- res$localdepth/nx
  result$depth <- res$depth/nx
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- max(result$depth)    
  result$num <- nx
  result$call <- match.call()
  result$tau <- tau
  result$use <- NA
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'halfspace'  
  class(result) <- 'localdepth'
  return(result)
}
