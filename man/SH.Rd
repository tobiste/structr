% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lund.R
\name{SH}
\alias{SH}
\title{Direction of maximum horizontal stress from the stress tensor}
\usage{
SH(S1, S2, S3, R, tol = .Machine$double.eps^0.5, ortho.tol = 1e-04)
}
\arguments{
\item{S1, S2, S3}{numeric three-column vectors. The principal stress orientations.
The variables hold the coordinates in the North, East and Down geographical
coordinate system, e.g. \code{S1 = c(s1N,s1E,s1D)}}

\item{R}{numeric. Relative magnitude of \code{S2} with respect to \code{S1} and \code{S3}:
\eqn{R = \frac{S1 - S2}{S1 - S3}}. Values ranging from 0 to 1, with 0 being
\code{S1==S2} and 1 being \code{S2==S3}.}

\item{tol}{Tolerance of comparison.}

\item{ortho.tol}{tolerance angle (in radians) for orthogonality check of the
three principal stress vectors.}
}
\value{
numeric angle in degrees. The direction of SH from North.
}
\description{
Calculates the direction of maximum horizontal stress using only the
directions of the principal stress and \eqn{R = \frac{S1 - S2}{S1 - S3}}.
This function Equations 11 and 10 from Lund and Townend (2007).
}
\examples{
# first example from https://www.snsn.se/SH/SHcode/benchmark.out
S1 <- Line(250.89, 70.07)
S3 <- Line(103.01, 17.07)
S2 <- vcross(S3, S1)
SH(S1, S2, S3, R = 1) #  70.89

R <- seq(0, 1, .05)
cbind(R, SH = sapply(R, function(x) {
  SH(S1, S2, S3, R = x)
}))
}
\references{
Lund and Townend, (2007). Calculating horizontal stress
orientations with full or partial knowledge of the tectonic stress tensor,
Geophys. J. Int., doi:\doi{10.1111/j.1365-246X.2007.03468.x}.
}
