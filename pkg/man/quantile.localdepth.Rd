\name{quantile.localdepth}
\alias{quantile.localdepth}
\title{Quantiles of the distribution of the size of simplices or ellipsoids}
\description{
The function measures the sizes of a sample from all the objects from a dataset according to a statistic and returns the corresponding quantiles. Depending on the method the objects are simplices or ellipsoids.
}
\usage{
quantile.localdepth(x, probs, use = c("volume", "diameter"), method = c("simplicial", "ellipsoid", "mahalanobis"), 
nsamp = "all", size = FALSE, dimension = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric; vector, dataframe or matrix}
  \item{probs}{numeric; vector of probabilities with values in [0,1]. Quantile orders of the statistic}
  \item{use}{character; the statistic used to measure the size of the objects. For \code{method} equal to \code{"simplicial"} or \code{"ellipsoid"} allowed statistics are \code{"volume"} and \code{"diameter"}. For \code{method} equal to \code{"mahalanobis"} this parameter is not used and the only available statistic is pairwise Mahalanobis' distance}
  \item{method}{character; the type of (local) depth to be evaluated}
  \item{nsamp}{character or numeric; the number of objects that are considered. If \code{"all"} then the size of all \code{choose(NROW(x), NCOL(x)+1)} objects is evaluated. Otherwise, a simple random sample with replacement of size \code{nsamp} is performed from the set of all possible objects. See details below}
  \item{size}{logical; If \code{TRUE} the size of all objects is also returned}
  \item{dimension}{numeric; only used with \code{method="ellipsoid"}. It is the squared length of the ellipsoid semimajor axis. If \code{dimension} is \code{NULL} then it is set to \code{NCOL(x)}}
  \item{\dots}{arguments passed to the next method}
}

\details{
When choosing between \code{nsamp='all'} or a numeric value for \code{nsamp} the value of \code{choose(NROW(x), NCOL(x)+1)} must be considered to avoid very long computing times. 
}

\value{
If \code{size} is \code{FALSE} then a vector with the quantiles is returned, otherwise a list with the following components:
  \item{quantile}{the quantiles}
  \item{stats}{the size of all objects}
  \item{call}{match.call()}
}
\author{Claudio Agostinelli and Mario Romanazzi}
\seealso{\code{\link{localdepth}}}
\examples{

set.seed(1234)
x <- matrix(rnorm(60, 0, 1), ncol=2)
volumesimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), size=TRUE)
volumesimplex$quantile
diametersimplex <- quantile.localdepth(x, probs=c(0.1, 0.2), size=TRUE, use='diameter')
diametersimplex$quantile
par(mfrow=c(2,2))
plot(ecdf(volumesimplex$stats), xlab='volume', main="ECDF of the simplices's volume")
boxplot(volumesimplex$stats, xlab='volume', main="Boxplot of the simplices's volume")
plot(ecdf(diametersimplex$stats), xlab='volume', main="ECDF of the simplices's diameter")
boxplot(diametersimplex$stats, xlab='volume', main="Boxplot of the simplices's diameter")

}

\keyword{multivariate}
\keyword{nonparametric}
