% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{classes}
\alias{classes}
\alias{Line}
\alias{Plane}
\alias{Fault}
\alias{as.line}
\alias{as.plane}
\alias{as.fault}
\alias{is.line}
\alias{is.plane}
\alias{is.fault}
\alias{is.spherical}
\title{Structure classes}
\usage{
Line(azimuth, plunge)

Plane(dip_direction, dip)

Fault(dip_direction, dip, azimuth, plunge, sense = NULL)

as.line(l)

as.plane(p)

as.fault(f)

is.line(l)

is.plane(p)

is.fault(f)

is.spherical(l)
}
\arguments{
\item{azimuth, plunge}{numeric vectors. Azimuth and plunge of a line (in
degrees)}

\item{dip_direction, dip}{numeric vectors. Dip direction and dip of a plane
(in degrees)}

\item{sense}{(optional) integer. Sense of the line on a fault plane. Either
\code{1}or \code{-1} for normal/dextral or thrust/sinistral offset, respectively.}

\item{l, p, f}{numeric vector or array containing the spherical coordinates
(1st element/column is azimuth, 2nd element/column is inclination, both in
degrees), or object of class \code{"line"}, \code{"plane"}, or \code{"fault"}}
}
\description{
\code{Line}, \code{Plane}, and \code{Fault} create a \code{"line"}, \code{"plane"}, and \code{"fault"}
class object, respectively, from the given set of values.

\code{as.line}, \code{as.plane}, and \code{as.fault} attempt to turn its argument into a
\code{"line"}, \code{"plane"}, and \code{"fault"} class object, respectively.

\code{is.line}, \code{is.plane}, and \code{is.fault} test if its argument is a \code{"line"},
\code{"plane"}, and \code{"fault"} class object, respectively.
}
\details{
\code{as.line}, \code{as.plane}, and \code{as.fault} return \code{TRUE} if \code{l}, \code{p}, and \code{f}
are an object of class \code{"line"}, \code{"plane"}, or \code{"fault"}, respectively, and
\code{FALSE} otherwise.

\code{is.spherical} returns \code{TRUE} if the argument's class is one of \code{"line"},
\code{"plane"}, or \code{"fault"} and \code{FALSE} otherwise

\code{as.line}, \code{as.plane}, and \code{as.fault} are is generic functions. If the
argument is a \code{"line"} or \code{"plane"} class, it will be converted.
}
\examples{
x <- Line(120, 50) # create line
is.line(x) # test if line
as.plane(x) # convert to plane
}
