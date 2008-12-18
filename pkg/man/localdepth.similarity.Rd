\name{localdepth.similarity}
\alias{localdepth.similarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
localdepth.similarity(x, y = NULL, tau, use = c("volume", "diameter"), method = c("simplicial", "ellipsoid", "mahalanobis"), type = c("exact", "approx"), nsamp = "all", nmax = 1, tol = 10^(-9), location = NULL, covariance = NULL, weight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{tau}{ ~~Describe \code{tau} here~~ }
  \item{use}{ ~~Describe \code{use} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{nsamp}{ ~~Describe \code{nsamp} here~~ }
  \item{nmax}{ ~~Describe \code{nmax} here~~ }
  \item{tol}{ ~~Describe \code{tol} here~~ }
  \item{location}{ ~~Describe \code{location} here~~ }
  \item{covariance}{ ~~Describe \code{covariance} here~~ }
  \item{weight}{ ~~Describe \code{weight} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y = NULL, tau, use = c("volume", "diameter"), method = c("simplicial", 
    "ellipsoid", "mahalanobis"), type = c("exact", "approx"), 
    nsamp = "all", nmax = 1, tol = 10^(-9), location = NULL, 
    covariance = NULL, weight = NULL) 
{
    use <- match.arg(use)
    method <- match.arg(method)
    type <- match.arg(type)
    if (method == "simplicial") {
        if (is.circular(x)) 
            localdepth.similarity.simp.circular(x = x, y = y, 
                tau = tau, use = use)
        if (type == "exact") 
            localdepth.similarity.simp(x = x, y = y, tau = tau, 
                use = use, weight = weight)
        else localdepth.similarity.simp.approx(x = x, y = y, 
            tau = tau, use = use, nsamp = nsamp, nmax = nmax, 
            tol = tol)
    }
    else if (method == "mahalanobis") {
        if (is.circular(x)) 
            stop("method 'mahalanobis' is not implemented for circular data")
        localdepth.similarity.mahalanobis(x = x, y = y, tau = tau, 
            nsamp = nsamp, nmax = nmax, location = location, 
            covariance = covariance, weight = weight)
    }
    else {
        if (is.circular(x)) 
            stop("method 'ellipsoid' is not implemented for circular data")
        localdepth.similarity.ellipsoid(x = x, y = y, tau = tau, 
            use = use, nsamp = nsamp, nmax = nmax, tol = tol)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
