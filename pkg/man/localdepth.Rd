\name{localdepth}
\alias{localdepth}
\title{Local Depth}
\description{
The function evaluates the depth and the local depth for a given set of points with respect the dataset. 
}
\usage{
localdepth(x, y = NULL, tau, use = c("volume", "diameter"), method = c("simplicial", "ellipsoid", "mahalanobis"), type = c("exact", "approx"), nsamp = "all", nmax = 1, tol=10^(-9), dimension=NULL)
}

\arguments{
  \item{x}{numeric. Vector, dataframe or matrix. If \code{x} is a \code{\link{circular}} vector than a circular version is used. Avoid ties in this dataset by little perturbation of the data before to use the function. The function only warnings for ties.}
  \item{y}{a vector, a dataframe or a matrix with the same number of columns as \code{x}, or \code{NULL}. If \code{NULL}, \code{x} is used}
  \item{tau}{the threshold value for the evaluation of the localdepth. Use function \code{\link{quantile.localdepth}} in order to fix \code{tau} using a quantile order of the dimension of the objects (depending of the method simplex, ellipsoid, etc.)}
  \item{use}{the statistic used to measure the size of the objects, so far the volume and the diameter are available}
  \item{method}{the type of (Local) Depth you want to evaluate}
  \item{type}{how to performance the evaluation of the membership. The \code{type="approx"} is available only for the \code{method="simplicial"}}
  \item{nsamp}{the number of objects that should be considered. If "all" then an exaustive scanning is performed}
  \item{nmax}{maximum fraction (in the range (0,1]) of objects to consider in case \code{nsamp} is not equal to \code{all}. If \code{1} means the same as the exaustive scanning (\code{choose(NROW(x), NCOL(x)+1)} for Simplicial and Ellipsoid Depth)}
  \item{tol}{it is a tolerance parameter that should be fixed depending the precision of your machine. This is used to ensure that the points on the border of the objects are considered inside it}
  \item{dimension}{this argument is used only when \code{method="ellipsoid"} and it is the square of the dimension of the ellipsoid to use. Defaulf is \code{NCOL(x)}}
}
\details{
  ~~ If necessary, more details than the description above ~~

}
\value{
  The function return an object of \code{\link{class}} \code{localdepth} with the following components:

  \item{localdepth}{the Local Depth for the given \code{tau}}
  \item{depth}{the Depth}
  \item{max.localdepth}{it is just \code{max(localdepth)}}
  \item{max.depth}{it is just \code{max(depth)}}
  \item{num}{vector with two elements. The first is the number of objects used for the evaluation of the Depth the second is the number of objects used for the evaluation of the Local Depth}
  \item{call}{match.call() result. Note that this is call from the internal function}
  \item{tau}{the threshold}
  \item{use}{the statistic used to measure the size of the objects}
  \item{tol}{the tolerance}
  \item{x}{as in input}
  \item{y}{as in input}
  \item{type}{as in input}
  \item{nsamp}{as in input}
  \item{method}{as in input}
}

\references{
C. Agostinelli and M. Romanazzi (2007) Local depth of univariate distributions, Working paper, Dipartimento di Statistica, Universita' Ca' Foscari, Venezia

C. Agostinelli and M. Romanazzi (2008)  Local depth of multidimensional data, Working paper n. 3/2008, Dipartimento di Statistica, Universita' Ca' Foscari, Venezia
}

\author{Claudio Agostinelli and Mario Romanazzi}

\seealso{\code{\link{quantile.localdepth}}, \code{\link{plot.localdepth}}}
\examples{





}

\keyword{multivariate}
\keyword{nonparametric}
