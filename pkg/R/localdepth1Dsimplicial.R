#############################################################
#
#	localdepth1Dsimplicial function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 14, 2007
#	Version: 0.2
#
#	Copyright (C) 2007 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth1Dsimplicial <- function(x, y, tau, use) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)

  if (use=='diameter') nuse <- 1
  if (use=='volume') nuse <- 2
  if (use=='spherical') nuse <- 3
  
  ## number of couples
  nc <- choose(nx, 2)
  result <- .Fortran("LD1DS",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.integer(nc),
    as.double(tau),
    as.integer(nuse),               
    localdepth=double(ny),               
    depth=double(ny),
    diameters=double(nc),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- NULL
  result[[1]] <- result[[2]] <- result[[3]] <- NULL
  result$areas <- result$diameters
  result$num <- nc
  return(result)
}

#############################################################
#
#	localdepth1Dtukey function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 8, 2007
#	Version: 0.1
#
#	Copyright (C) 2007 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

####### DA FARE

localdepth1Dtukey <- function(x, y=NULL, use=c('R', 'fortran')) {
  use <- match.arg(use)
  nx <- length(x)
  if (is.null(y)) {
     y <- x
  }
  ny <- length(y)
  differenze <- rep(x,times=length(y))-rep(y,each=length(x))
  diametri <- unique(abs(differenze))
  nc <- length(diametri)
  if (use=='fortran') {
    result <- .Fortran("LD1DT",
	as.double(x),
	as.double(y),
	as.integer(nx), 
	as.integer(ny),
	as.integer(nc),
	depth=mat.or.vec(nc,ny),
	diameters=double(nc),
	PACKAGE = "localdepth")
  } else {
    depth <- matrix(0, nrow=nc, ncol=ny)
    for (i1 in 1:ny) {
      temp1 <- differenze[(nx*(i1-1)+1):(nx*i1)]
      for (i2 in 1:nc) {
          temp2 <- temp1[abs(temp1)<=diametri[i2]]
          depth[i2, i1] <- min(sum(temp2>=0), sum(temp2<=0))
      }
    }
    result <- list()
    result$depth <- depth

  }
  result$diameters <- diametri 
  result$areas <- result$diameters
  return(result)
}
