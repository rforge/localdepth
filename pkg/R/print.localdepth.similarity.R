#############################################################
#
#	print.localdepth.similarity function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 2, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
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
    cat("Threshold parameter 'tau' is set to: ", x$tau, "\n")
    cat("The dimension is measured by: ", x$use, "\n")
    cat("Membership evaluation: ", x$type, "\n")
    if (x$nsamp=='all')
      cat("All objects were explored \n")
    else
      cat("A Monte Carlo approach was used on: ", x$nun[1], " objects in order to have ", x$nun[1] ,"\n")
    cat("\n")
    invisible(x)
}
