#############################################################
#
#	print.localdepth.similarity function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: April, 17, 2015
#	Version: 0.2-1
#
#       Copyright (C) 2015 Claudio Agostinelli
#	Copyright (C) 2008-2014 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

print.localdepth.similarity <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Local Depth Similarity Statistics: \n")
    print(summary(x$localdepth[upper.tri(x$localdepth, diag=TRUE)]), digits=digits, ...)
    cat("Depth Similarity Statistics: \n")
    print(summary(x$depth[upper.tri(x$depth, diag=TRUE)]), digits=digits, ...)
    cat("\n")
    cat("Information \n")
    cat("Method: ", x$method, "\n")
    if (is.function(x$tau))
      cat("Threshold parameter 'tau' is a function\n")    
##      cat("Function 'tau' is set to: ", substitute(x$tau), "\n")
    else if (length(tau)==1)
      cat("Threshold parameter 'tau' is set to: ", x$tau, "\n")
    else
      cat("Threshold parameter 'tau' is a vector, a summary is: \n", summary(x$tau), "\n")
    if (!is.na(x$use))
      cat("The dimension is measured by: ", x$use, "\n")
    if (!is.na(x$type))    
      cat("Membership evaluation: ", x$type, "\n")    
    if (x$nsamp=='all')
      cat("All objects were explored \n")
    else {
      if (x$method!='mahalanobis') {
        cat("A Monte Carlo approach was used on: ", x$num[1], " objects in order to have ", x$num[2])
        if (x$method=='simplicial')
          cat(" simplicies ")
        else if (x$method=='ellipsoid')
          cat(" ellipsoids ")
        cat("with size smaller than 'tau'\n")
      } else {
        cat("A Monte Carlo approach was used on: ", x$num[1], " objects\n")
      }
    }
    cat("\n")
    invisible(x)
}
