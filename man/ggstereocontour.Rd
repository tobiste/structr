% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_stereonet.R
\name{ggstereocontour}
\alias{ggstereocontour}
\alias{geom_contour_stereo}
\alias{geom_contourf_stereo}
\title{Stereonet contouring using ggplot}
\usage{
geom_contour_stereo(
  data,
  ngrid = 200,
  hw = NULL,
  optimal_bw = c("cross", "rot"),
  norm = FALSE,
  threshold = 0,
  ...
)

geom_contourf_stereo(
  data,
  ngrid = 200,
  hw = NULL,
  optimal_bw = c("cross", "rot"),
  norm = FALSE,
  smooth = FALSE,
  threshold = 0,
  ...
)
}
\arguments{
\item{data}{data.frame containing the orientation}

\item{ngrid}{integer. Resolution of density calculation.}

\item{hw}{numeric. Kernel bandwidth in degree.}

\item{optimal_bw}{character. Calculates an optimal kernel bandwidth
using the cross-validation algorithm (\code{'cross'}) or the rule-of-thumb (\code{'rot'})
suggested by Garcia-Portugues (2013). Ignored when \code{hw} is specified.}

\item{norm}{logical. Should the densities be normalized?}

\item{threshold}{numeric. Cut-off for low densities.}

\item{...}{arguments passed to \code{\link[ggplot2:geom_contour]{ggplot2::geom_contour()}}, \code{\link[ggplot2:geom_contour]{ggplot2::geom_contour_filled()}}, or \code{\link[ggplot2:geom_tile]{ggplot2::geom_tile()}}}

\item{smooth}{logical. Whether \code{\link[ggplot2:geom_tile]{ggplot2::geom_tile()}} should be used for plotting.}
}
\value{
ggplot
}
\description{
Stereonet contouring using ggplot
}
\examples{
if (require("mapproj")) {
  test_data <- rbind(
    rvmf(100, mu = Line(90, 45), k = 10),
    rvmf(50, mu = Line(0, 0), k = 20)
  ) |> as.line()

  ggstereo() +
    geom_contourf_stereo(gg(test_data)) +
    ggplot2::scale_fill_viridis_d(option = "A") +
    # guides(fill = guide_colorsteps(barheight = unit(8, "cm"), show.limits = TRUE)) +
    geom_contour_stereo(gg(test_data), color = "grey") +
    ggplot2::geom_point(data = gg(test_data), ggplot2::aes(x = x, y = y), color = "lightgrey") +
    ggframe()

  ggstereo() +
    geom_contourf_stereo(gg(test_data), norm = TRUE, bins = 50, threshold = .1) +
    ggplot2::scale_fill_viridis_d(option = "A")
}
}
\references{
Garcia Portugues, E. (2013). Exact risk improvement of
bandwidth selectors for kernel density estimation with directional data.
Electronic Journal of Statistics, 7, 1655-1685.
}
