% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sph_distribution.R
\name{dist.mle}
\alias{dist.mle}
\alias{kent.mle}
\alias{vmf.mle}
\title{MLE of spherical rotational symmetric distributions}
\source{
Adapted from \code{\link[Directional:kent.mle]{Directional::kent.mle()}} and \code{\link[Directional:vmf.mle]{Directional::vmf.mle()}}
}
\usage{
kent.mle(x)

vmf.mle(x)
}
\arguments{
\item{x}{numeric. Can be three element vector, three column array, or an
object of class \code{"line"} or \code{"plane"}}
}
\description{
Estimates the parameters of a von Mises-Fisher or Kent distribution.
}
\examples{
x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
kent.mle(x)
vmf.mle(x)
}
