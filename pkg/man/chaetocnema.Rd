\name{chaetocnema}
\alias{chaetocnema}
\docType{data}
\title{Chaectonema dataset}
\description{
The Chaectonema dataset  are anatomical measures of 74 specimens of fleabeetles, of three different genera, Ch. concinna (A, 21 cases), Ch. heptapotamica (B, 31 cases) and Ch. heikertingeri (C, 22 cases).
}
\usage{data(chaetocnema)}
\format{
  A data frame with 74 observations on the following 3 variables.
  \describe{
    \item{\code{V1}}{}
    \item{\code{V2}}{}
    \item{\code{lab}}{a factor with levels \code{A}, \code{B}, \code{C}, see the description above}
  }
}
\details{

}
\source{
A.A. Lubischew (1962) On the use of discriminant functions in taxonomy. Biometrics, 18:4, 455--477.
}
\references{
  ~~ possibly secondary sources and usages ~~
}

\examples{
data(chaetocnema)
lab <- chaetocnema[,3]
chaetocnema <- chaetocnema[,-3]
set.seed(1234)
chaetosd <- apply(chaetocnema, 2, sd)
chaetocnema <- chaetocnema+cbind(rnorm(nrow(chaetocnema), 0, chaetosd[1]/100), rnorm(nrow(chaetocnema), 0, chaetosd[2]/100))


gx <- seq(min(chaetocnema[,1]), max(chaetocnema[,1]), length=60)
gy <- seq(min(chaetocnema[,2]), max(chaetocnema[,2]), length=60)
y <- as.matrix(expand.grid(gx, gy))

tauvolume <- quantile.localdepth(chaetocnema, probs=0.20, use='volume')
system.time(depthcontours <- localdepth(x=chaetocnema, y=y, tau=tauvolume, use='volume', method='simplicial'))
system.time(depthcontoursapprox <- localdepth(x=chaetocnema, y=y, tau=tauvolume, use='volume', method='simplicial', type='approx'))

par(mfrow=c(2,2))
plot(chaetocnema, pch=as.character(lab), main='Exact Simplicial Depth')
contour(x=gx, y=gy, z=matrix(depthcontours$depth/depthcontours$max.depth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Approx Simplicial Depth')
contour(x=gx, y=gy, z=matrix(depthcontoursapprox$depth/depthcontoursapprox$max.depth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Exact Simplicial Local Depth, Tau=20%')
contour(x=gx, y=gy, z=matrix(depthcontours$localdepth/depthcontours$max.localdepth, nrow=length(gx), ncol=length(gy)), add=TRUE)

plot(chaetocnema, pch=as.character(lab), main='Approx Simplicial Local Depth, Tau=20%')
contour(x=gx, y=gy, z=matrix(depthcontoursapprox$localdepth/depthcontoursapprox$max.localdepth, nrow=length(gx), ncol=length(gy)), add=TRUE)
}
\keyword{datasets}
