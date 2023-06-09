% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sph_distribution.R
\name{rfb}
\alias{rfb}
\title{Spherical Fisher-Bingham distribution}
\source{
Adapted from \code{\link[Directional:rfb]{Directional::rfb()}}
}
\usage{
rfb(n = 100, mu = c(1, 0, 0), k = 5, A)
}
\arguments{
\item{n}{integer. number of random samples to be generated}

\item{mu}{numeric. Can be three element vector, three column array, or an
object of class \code{"line"} or \code{"plane"}}

\item{k}{numeric. The concentration parameter (\eqn{\kappa})}

\item{A}{symmetric matrix}
}
\description{
Simulation of random values from a spherical Fisher-Bingham distribution.
}
\examples{
\dontrun{
x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
stereoplot()
stereo_point(x)
}
}
