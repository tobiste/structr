% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faults.R
\name{fault_analysis}
\alias{fault_analysis}
\title{Simple fault analysis}
\usage{
fault_analysis(x, ptangle = 90)
}
\arguments{
\item{x}{fault object}

\item{ptangle}{angle between P and T axes in degrees (90 by default) .}
}
\value{
list
}
\description{
Calculates PT-axes, kinematic plane (M), and dihedra separation plane (d)
}
\examples{
f <- Fault(120, 60, 110, 58, 1)
fault_analysis(f)
}
