\name{quantile.localdepth}
\alias{quantile.localdepth}
\title{Quantiles of the distribution of the size of objects}
\description{
The function measures all (or resampling from all) the objects constructed from a dataset according to a statistic and return the requested quantiles of the distribution
}
\usage{
quantile.localdepth(x, probs, use = c("volume", "diameter"), method = c("simplicial", "ellipsoid", "mahalanobis"), 
nsamp = "all", all = FALSE, dimension = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric. Vector, dataframe or matrix}
  \item{probs}{numeri.: Vector of probabilities with values in [0,1]. These are the order of the quantiles requested}
  \item{use}{character. The statistic used to measures the size of the objects}
  \item{method}{character. The type of (Local) Depth you want to evaluate}
  \item{nsamp}{character or numeric: the number of objects that should be considered. If "all" then an exaustive scanning is performed}
  \item{all}{logical. If \code{TRUE} the size of all objects is also returned}
  \item{dimension}{numeric. This argument is used only when \code{method="ellipsoid"} and it is the square of the dimension of the ellipsoid. If \code{dimension} is \code{NULL} then it is set to \code{NCOL(x)}}
  \item{\dots}{arguents passed to the next method}
}

\details{


}

\value{
If \code{all} is \code{FALSE} then a vector with the quantiles is returned otherwise a list with the follows components is returned:
  \item{quantile}{the quantiles}
  \item{stats}{the size of all objects}
  \item{call}{match.call()}
}
\author{Claudio Agostinelli and Mario Romanazzi}
\seealso{\code{\link{localdepth}}}
\examples{

set.seed(1234)
x <- matrix(rnorm(60, 0, 1), ncol=2)
volumesimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), all=TRUE)
volumesimplex$quantile
diametersimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), all=TRUE, use='diameter')
diametersimplex$quantile
par(mfrow=c(2,2))
plot(ecdf(volumesimplex$stats), xlab='volume', main="ECDF of the simplices's volume")
boxplot(volumesimplex$stats, xlab='volume', main="Boxplot of the simplices's volume")
plot(ecdf(diametersimplex$stats), xlab='volume', main="ECDF of the simplices's diameter")
boxplot(diametersimplex$stats, xlab='volume', main="Boxplot of the simplices's diameter")

}

\keyword{multivariate}
\keyword{nonparametric}
