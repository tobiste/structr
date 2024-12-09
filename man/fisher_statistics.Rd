% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/math.R
\name{fisher_statistics}
\alias{fisher_statistics}
\title{Fisher's statistics}
\usage{
fisher_statistics(x, w = NULL, p = 0.05)
}
\arguments{
\item{x}{numeric. Can be three element vector, three column array, or an
object of class \code{"line"} or \code{"plane"}}

\item{w}{numeric. Weights}

\item{p}{numeric. Significance level (\code{0.05} by default), corresponding to
\eqn{100 * (1-p) - 95\%} confidence level.}
}
\value{
list, with
\describe{
\item{\code{"k"}}{estimated concentration parameter \eqn{\kappa} for the von Mises-Fisher
distribution}
\item{\code{"csd"}}{estimated angular standard deviation enclosing 63\% of the orientation data. Angle is in degrees if \code{x} is a spherical object, and raidan if otherwise.}
\item{\code{"a95"}}{Confidence limit for given \code{p}. Angle is in degrees if \code{x} is a spherical object, and raidan if otherwise.}
}
}
\description{
Estimates concentration parameter, angular standard deviation, and
confidence limit.
}
\examples{
x <- rvmf(100, mu = Line(120, 50), k = 5)
fisher_statistics(x)
}
