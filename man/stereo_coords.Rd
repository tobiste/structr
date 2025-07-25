% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stereonet.R
\name{stereo_coords}
\alias{stereo_coords}
\title{Stereographic projection}
\usage{
stereo_coords(az, inc, upper.hem = FALSE, earea = TRUE, r = 1)
}
\arguments{
\item{az, inc}{numeric vectors. Azimuth and Inclination in degrees.}

\item{upper.hem}{logical. Whether the projection is shown for upper
hemisphere (\code{TRUE}) or lower hemisphere (\code{FALSE}, the default).}

\item{earea}{logical \code{TRUE} for Lambert equal-area projection (also "Schmidt net"; the default), or
\code{FALSE} for meridional stereographic projection (also "Wulff net" or "Stereonet").}

\item{r}{numeric. Radius of circle. Default is \code{1} for unit circle.}
}
\value{
two-column vector with the transformed coordinates
}
\description{
Transformation of spherical coordinates into the stereographic projection
}
\examples{
stereo_coords(90, 10)
stereo_coords(90, 10, earea = TRUE, upper.hem = TRUE)
}
