# library(stereonet)
# library(pracma)

# functions are from mplstereonet

# sph2cart <- function(lon, lat) {
#   x <- cos(lat) * cos(lon)
#   y <- sin(lat) * sin(lon)
#   z <- sin(lat)
#   cbind(x, y, z)
# }
#
# cart2sph <- function(vec) {
#   x <- vec[, 1]
#   y <- vec[, 2]
#   z <- vec[, 3]
#   r <- sqrt(x**2 + y**2 + z**2)
#   lat <- asin(z / r)
#   lon <- atan2(y, x)
#   cbind(lon, lat)
# }

#' Point density
#'
#' This function actually calculates the point density of the input
#' (`lons` and `lats`) points at a series of counter stations. Creates
#' `gridsize` regular grid of counter stations in lat-long space, calculates the
#'  distance to all input points at each counter station, and then calculates
#'  the density using `FUN`.  Each input point is weighted by the corresponding
#'  item of `weights`. The weights are normalized to 1 before calculation.
#'
#' @param azi,inc degree
#' @param FUN The method of density estimation to use. Defaults to [exponential_kamb()].
#' @param sigma (optional) numeric. The number of standard deviations defining the expected number of
#' standard deviations by which a random sample from a uniform
#' distribution of points would be expected to vary from being evenly
#' distributed across the hemisphere.  This controls the size of the
#' counting circle, and therefore the degree of smoothing.  Higher `sigma`s
#' will lead to more smoothing of the resulting density distribution. This
#' parameter only applies to Kamb-based methods.  Defaults to `3`.
#' @param n (optional) two element vector. The size of the grid that the density is estimated on.
#' @param weights (optional) numeric vector of length of `lons` and `lats`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#'
#' @returns list
#' @export
count_points <- function(azi, inc, FUN = exponential_kamb, sigma = 3, n = 100, weights, ...) {
  # Ensure lons and lats are numeric vectors
  # lons <- as.numeric(lons)
  # lats <- as.numeric(lats)

  # Grid setup
  # bound <- pi / 2
  # nrows <- gridsize[1]
  # ncols <- gridsize[2]
  # # grid_lon <- matrix(seq(-bound, bound, length.out = ncols), nrow = ncols, ncol = nrows)
  # # grid_lat <- matrix(seq(-bound, bound, length.out = nrows), nrow = ncols, ncol = nrows)
  # grid_azi <- seq(0, 2*pi, length.out = ncols)
  # grid_inc <- seq(0, pi, length.out = nrows)
  # grid <- expand.grid(azi = grid_azi, inc = grid_inc) |> rad2deg()

  grid <- blank_grid(n = n, ...)

  # Stereonet math transformations to cartesian coordinates
  # xyz_counters <- lin2vec0(grid[, 1], grid[, 2])
  xyz_points <- lin2vec0(azi, inc)

  # totals <- numeric(nrow(grid$grid))

  for (i in seq_along(grid$grid[, 1])) {
    cos_dist <- abs(grid$grid[i, ] %*% t(xyz_points))
    density_scale <- do.call(FUN, args = list(cos_dist, sigma))
    density <- density_scale$count * weights
    scale <- density_scale$units
    # totals[i] <- (sum(density) - 0.5) / scale
    grid$density[i] <- (sum(density) - 0.5) / scale
  }

  # totals[totals < 0] <- 0
  grid$density[grid$density < 0] <- 0

  # xyz_counters_sph <- cart2sph(xyz_counters)
  # counter_lon <- xyz_counters_sph[, 1]
  # counter_lat <- xyz_counters_sph[, 2]

  # list(
  #   xi = matrix(counter_lon, gridsize),
  #   yi = matrix(counter_lat, gridsize),
  #   #lon = grid_lon,
  #   #lat = grid_lat,
  #   density = matrix(totals, gridsize)
  #   #zi = totals
  # )
  # cbind(xyz_counters, density=totals) # return as Cartesian
  # cbind(grid, density = totals) # return in lon-lat space (degrees)
  grid
}

#' Estimates point density
#'
#' Estimates point density of the given linear orientation measurements
#' Returns a regular (in lat-long space) grid of density estimates over a
#' hemispherical surface.
#'
#' @param x spherical object
#' @param n (optional) numeric. The size of the grid that the density is estimated on. Defaults to `100`.
#' @param arguments passed to [count_points()]
#'
#'
#' @returns list with the cartesian coordinates of the grid and density values of the regularly gridded density estimates.
#' @export
#'
#' @examples
#' x <- Line(c(120, 315, 86), c(22, 85, 31))
#' density_grid(x, n = 100)
density_grid <- function(x, FUN = exponential_kamb, sigma = 3, n = 100, weights = NULL, ...) {
  # if (length(gridsize) == 1) gridsize <- c(gridsize, gridsize)

  if (!is.spherical(x)) x <- to_spherical(x)

  if (!is.line(x)) x <- as.line(x)

  azi <- x[, 1]
  inc <- x[, 2]

  if (is.null(weights)) weights <- rep(1, nrow(x))

  # normalize weights to 1
  # weights <- as.numeric(weights) / mean(weights)
  weights <- as.numeric(weights) / max(weights)

  count_points(azi, inc, FUN = FUN, sigma = sigma, n = n, weights = weights, ...)
}


kamb_radius <- function(n, sigma) {
  a <- sigma^2 / (n + sigma^2)
  return(1 - a)
}

kamb_units <- function(n, radius) {
  return(sqrt(n * radius * (1 - radius)))
}


#' Density estimation
#'
#' The methods of density estimation in a stereonet.
#'
#' @param cos_dist cosine distances
#' @param sigma (optional) numeric. The number of standard deviations defining the expected number of
#' standard deviations by which a random sample from a uniform
#' distribution of points would be expected to vary from being evenly
#' distributed across the hemisphere.  This controls the size of the
#' counting circle, and therefore the degree of smoothing.  Higher sigmas
#' will lead to more smoothing of the resulting density distribution. This
#' parameter only applies to Kamb-based methods.  Defaults to 3.
#'
#' @details [exponential_kamb()]: Kamb with exponential smoothing
#' A modified Kamb method using exponential smoothing (ref1). Units are
#' in numbers of standard deviations by which the density estimate
#' differs from uniform.
#'
#' [linear_kamb()]: Kamb with linear smoothing
#' A modified Kamb method using linear smoothing (ref1).  Units are in
#' numbers of standard deviations by which the density estimate differs from uniform.
#'
#' [kamb()]: Kamb with no smoothing
#' Kamb's method (ref2) with no smoothing. Units are in numbers of standard
#' deviations by which the density estimate differs from uniform.
#'
#' [schmidt()]: 1% counts. The traditional "Schmidt" (a.k.a. 1%) method.
#' Counts points within a counting circle comprising 1% of the total area of the
#' hemisphere. Does not take into account sample size.  Units are in points per 1% area.
#'
#' @references
#' Vollmer, 1995. C Program for Automatic Contouring of Spherical
#' Orientation Data Using a Modified Kamb Method. Computers &
#'   Geosciences, Vol. 21, No. 1, pp. 31--49.
#'
#'   Kamb, 1959. Ice Petrofabric Observations from Blue Glacier,
#'   Washington, in Relation to Theory and Experiment. Journal of
#'   Geophysical Research, Vol. 64, No. 11, pp. 1891--1909.
#'
#' @returns list
#' @name density-funs
NULL

#' @rdname density-funs
#' @export
exponential_kamb <- function(cos_dist, sigma = 3) {
  n <- length(cos_dist)
  f <- 2 * (1 + n / sigma^2)
  count <- exp(f * (cos_dist - 1))
  units <- sqrt(n * (f / 2 - 1) / f^2)
  return(list(count = count, units = units))
}

#' @rdname density-funs
#' @export
linear_inverse_kamb <- function(cos_dist, sigma = 3) {
  n <- length(cos_dist)
  radius <- kamb_radius(n, sigma)
  f <- 2 / (1 - radius)
  cos_dist <- cos_dist[cos_dist >= radius]
  count <- f * (cos_dist - radius)
  return(list(count = count, units = kamb_units(n, radius)))
}

#' @rdname density-funs
#' @export
square_inverse_kamb <- function(cos_dist, sigma = 3) {
  n <- length(cos_dist)
  radius <- kamb_radius(n, sigma)
  f <- 3 / (1 - radius)^2
  cos_dist <- cos_dist[cos_dist >= radius]
  count <- f * (cos_dist - radius)^2
  return(list(count = count, units = kamb_units(n, radius)))
}

#' @rdname density-funs
#' @export
kamb_count <- function(cos_dist, sigma = 3) {
  n <- length(cos_dist)
  dist <- kamb_radius(n, sigma)
  count <- ifelse(cos_dist >= dist, 1, 0)
  return(list(count = count, units = kamb_units(n, dist)))
}

#' @rdname density-funs
#' @export
schmidt_count <- function(cos_dist, sigma = NULL) {
  radius <- 0.01
  count <- ifelse((1 - cos_dist) <= radius, 1, 0)
  count <- 0.5 / length(count) + count
  return(list(count = count, units = length(cos_dist) * radius))
}

reshape_grid <- function(m, n) {

}

#' @importFrom graphics contour filled.contour
#' @importFrom stats xtabs
stereo_density <- function(x, nlevels = 20, ..., filled = FALSE, upper.hem = FALSE) {
  d <- density_grid(x, ...)
  d$grid <- fix_symm(d$grid)

  d_sph <- vec2line(d$grid) |> cbind(density = d$density)

  m <- xtabs(density ~ azimuth + plunge, data = d_sph)

  if (filled) {
    filled.contour(as.numeric(rownames(m)), as.numeric(colnames(m)), m,
      nlevels = nlevels, color.palette = viridis::magma,
      # plot.axes = {
      #   axis(1); axis(2)
      #
      #   coords = expand.grid(x, y)
      #   projected = stereo_coords(coords[, 1], coords[, 2], upper.hem = upper.hem)
      #   points(projected[, 1], projected[, 2], col = viridis::magma(10)[cut(z, 10)])
      # }
    )
  } else {
    contour(m, nlevels = nlevels)
  }
}

fix_symm <- function(x) {
  x[x[, 3] < 0, ] <- v_antipode(x[x[, 3] < 0, ])
  x
}



#' Plot density grid in a stereonet
#'
#' @param x spherical data
#' @param ... arguments passed top [blank_grid()]
#' @param pal color function
#' @param binned logical. Whether the density colors should be binned
#' @param breaks for binning
#'
#' @return plot
#' @export
#'
#' @examples
#' test <- Line(c(90, 85, 105), c(45, 50, 40))
#' stereoplot(guides = F, centercross = FALSE)
#' stereo_density_grid(test, binned = FALSE, n = 10000)
#' stereo_point(test, cex = 2, col = "red")
stereo_density_grid <- function(x, ..., binned = FALSE, pal = viridis::viridis, breaks = 5) {
  d <- density_grid(x, ...)
  d$grid <- fix_symm(d$grid)

  if (binned) {
    cols <- bin_color(d$density, breaks = breaks, pal = pal)
  } else {
    cols <- assign_col(d$density, pal = pal)
  }

  stereo_point(to_spherical(d$grid), col = cols)
}
