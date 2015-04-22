#############################################################
#
#	depth2Dsimplicialsimilarity function
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: February, 20, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#
#############################################################

depth2Dsimplicialsimilarity <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  result <- list()

  ## numero di terne
  nt <- choose(nx, 3)
  depth <- rep(0, ny*nt)
  res <- .C("d2Dsimpsim",
            x1 = as.double(x[,1]),
            x2 = as.double(x[,2]),
            y1 = as.double(y[,1]),
            y2 = as.double(y[,2]),
            nx = as.integer(nx),
            ny = as.integer(ny),
            nt = as.integer(nt),
            depth = as.double(depth),
            DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")
  res$depth <- matrix(res$depth, nrow=ny, ncol=nt, byrow=FALSE)
  res$depth <- res$depth%*%t(res$depth)
  result$depth <- res$depth
  result$call <- match.call()
  result$num <- nt
  class(result) <- 'localdepth.similarity'
  return(result)
}
