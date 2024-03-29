% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mohR.R
\name{PR_stress}
\alias{PR_stress}
\title{Principal stresses from 2D stress components}
\usage{
PR_stress(sx, sy, txy)
}
\arguments{
\item{sx}{numeric. Normal stress acting on plane facing in X direction
(\deqn{sigma_x}).}

\item{sy}{numeric. Normal stress acting on plane facing in Y direction
(\deqn{\sigma_y}).}

\item{txy}{numeric. Shear stress acting on planes facing X and Y
(\deqn{\tau_{xy}}).}
}
\value{
angle in degrees
}
\description{
Determines the principal stresses and their orientations from the stress
components \deqn{sigma_x}, \deqn{\sigma_y}, \deqn{\tau_{xy}}.
}
\references{
Richard J. Lisle (1999)
}
