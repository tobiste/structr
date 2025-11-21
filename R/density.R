#' @keywords internal
blank_grid_regular <- function(n, r = 1) {
  x_grid <- seq(-1, 1, length.out = n)
  y_grid <- seq(-1, 1, length.out = n)
  grid_schmidt <- expand.grid(x_grid, y_grid) |> as.matrix()

  grid_cart <- .schmidt2cart(grid_schmidt[, 1], grid_schmidt[, 2])
  values <- rep(0, times = nrow(grid_cart))
  list(grid = grid_cart, density = values)
}

#' @keywords internal
.schmdit2spherical <- function(x, y, r = 1) {
  tq <- sqrt(x^2 + y^2)

  # Recover inclination in radians
  inc_rad <- 2 * (pi / 4 - asin(tq / (r * sqrt(2))))

  # Recover azimuth in radians
  az_rad <- atan2(x, y) %% (2 * pi)

  # Return as matrix with named columns
  cbind(az_rad = az_rad, inc_rad = inc_rad)
}

#' @keywords internal
.schmidt2cart <- function(x, y, r = 1) {
  res <- .schmdit2spherical(x, y, r) |> rad2deg()
  Line(res[, 1], res[, 2]) |>
    Vec3() |>
    unclass()
}


#' Point density
#'
#' This function calculates the point density of the input
#' spherical points at a series of counter stations. Creates
#' `gridsize` regular grid of counter stations, calculates the
#'  distance to all input points at each counter station, and then calculates
#'  the density using `FUN`.  Each input point is weighted by the corresponding
#'  item of `weights`. The weights are normalized to 1 before calculation.
#'
#' @param azi,inc degrees.
#' @param FUN The method of density estimation to use. Defaults to [exponential_kamb()].
#' @param sigma (optional) numeric. The number of standard deviations defining the expected number of
#' standard deviations by which a random sample from a uniform
#' distribution of points would be expected to vary from being evenly
#' distributed across the hemisphere.  This controls the size of the
#' counting circle, and therefore the degree of smoothing.  Higher `sigma`s
#' will lead to more smoothing of the resulting density distribution. This
#' parameter only applies to Kamb-based methods.  Defaults to `3`.
#' @param ngrid numeric. The size of the grid that the density is estimated on.
#' @param weights (optional) numeric vector of length of `azi`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#' @param r numeric. radius of stereonet circle
#'
#' @returns list
count_points <- function(azi, inc, FUN, sigma, ngrid, weights, r) {
  # Grid setup
  grid <- blank_grid_regular(n = ngrid, r = r)
  grid_coords <- grid$grid

  # Stereonet math transformations to Cartesian coordinates
  xyz_points <- Line(azi, inc) |>
    Vec3() |>
    unclass()

  for (i in seq_along(grid_coords[, 1])) {
    cos_dist <- abs(grid_coords[i, ] %*% t(xyz_points))
    density_scale <- FUN(cos_dist, sigma)
    density <- density_scale$count * weights
    scale <- density_scale$units
    grid$density[i] <- (sum(density) - 0.5) / scale
  }

  grid$density[grid$density < 0] <- 0

  grid
}


#' Estimates point density
#'
#' Estimates point density of the given linear orientation measurements
#' Returns a regular (in lat-long space) grid of density estimates over a
#' hemispherical surface.
#'
#' @inheritParams sph_mean
#' @param weights (optional) numeric vector of length of `azi`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param kamb logical. Whether to use the von Mises-Fisher kernel density estimation (`FALSE`) or Kamb's method (`TRUE`, the default).
#' @param ... arguments passed to [count_points()]
#'
#'
#' @returns list with the cartesian coordinates of the grid and density values of the regularly gridded density estimates.
#' @noRd
#'
#' @examples
#' x <- rvmf(100, mu = Line(120, 50))
#' res <- density_grid(x, FUN = kamb_count, ngrid = 100, sigma = 4)
#' lapply(res, head)
#'
#' res2 <- density_grid(x, kamb = FALSE, kernel_method = "cross")
#' lapply(res2, head)
density_grid <- function(x, weights = NULL, upper.hem = FALSE, kamb = TRUE, ...) {
  if (!is.Line(x)) x <- Line(x)

  azi <- x[, 1]
  if (upper.hem) {
    azi <- azi + 180
  }
  inc <- x[, 2]

  if (is.null(weights)) weights <- rep(1, nrow(x))

  # normalize weights to 1
  # weights <- as.numeric(weights) / mean(weights)
  weights <- as.numeric(weights) / max(weights)

  if (kamb) {
    count_points(azi, inc, weights = weights, ...)
  } else {
    res <- vmf_kerncontour(.full_hem(azi, inc), ...)
    # grid <- expand.grid(Lat = res$lat - 90, Long = res$long - 180)
    grid <- expand.grid(inc = res$lat - 90, azi = res$long - 180)

    list(
      grid = Line(grid$azi, grid$inc) |> Vec3(),
      density = c(res$den)
    )
  }
}

#' @keywords internal
kamb_radius <- function(n, sigma) {
  a <- sigma^2 / (n + sigma^2)
  return(1 - a)
}

#' @keywords internal
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
#' [linear_inverse_kamb()]: Kamb with linear smoothing
#' A modified Kamb method using linear smoothing (ref1).  Units are in
#' numbers of standard deviations by which the density estimate differs from uniform.
#'
#' [square_inverse_kamb()]: Kamb with squared smoothing
#' A modified Kamb method using squared smoothing (ref1).  Units are in
#' numbers of standard deviations by which the density estimate differs from uniform.
#'
#' [kamb_count()]: Kamb with no smoothing
#' Kamb's method (ref2) with no smoothing. Units are in numbers of standard
#' deviations by which the density estimate differs from uniform.
#'
#' [schmidt_count()]: 1% counts. The traditional "Schmidt" (a.k.a. 1%) method.
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
#' @importFrom Directional euclid vmfkde.tune vmf.mle
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
# not working
linear_inverse_kamb <- function(cos_dist, sigma = 3) {
  n <- length(cos_dist)
  radius <- kamb_radius(n, sigma)
  f <- 2 / (1 - radius)
  cos_dist <- cos_dist[cos_dist >= radius]
  count <- f * (cos_dist - radius)
  return(list(count = count, units = kamb_units(n, radius)))
}

#' @rdname density-funs
# not working
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
  count <- 1 / (2 * length(count)) + count
  return(list(count = count, units = length(cos_dist) * radius))
}

#' @keywords internal
vmf_kerncontour <- function(u, hw = NULL, kernel_method = c("cross", "rot"), ngrid = 100) {
  n <- nrow(u)
  x <- Directional::euclid(u)

  if (is.null(hw)) {
    kernel_method <- match.arg(kernel_method, c("cross", "rot"))

    if (kernel_method == "cross") {
      hw <- as.numeric(Directional::vmfkde.tune(x, low = 0.1, up = 1)[1])
    } else {
      k <- Directional::vmf.mle(x, fast = TRUE)$kappa
      hw <- ((8 * sinh(k)^2) / (k * n * ((1 + 4 * k^2) * sinh(2 * k) - 2 * k * cosh(2 * k))))^(1 / 6)
    }
  } else {
    hw <- deg2rad(hw)
  }

  kappa_val <- 1 / (hw^2)
  # cpk <- 1 / (sqrt(hw^2) * (2 * pi)^1.5 * besselI(kappa_val, 0.5))
  cpk <- 1 / (sqrt(hw^2) * (2 * pi)^1.5 * besselI(kappa_val, 0.5, expon.scaled = TRUE) * exp(kappa_val))

  lat_grid <- seq(0, 180, length.out = ngrid)
  long_grid <- seq(0, 360, length.out = ngrid)
  den_mat <- matrix(nrow = ngrid, ncol = ngrid)
  for (i in 1:ngrid) {
    for (j in 1:ngrid) {
      y <- Directional::euclid(c(lat_grid[i], long_grid[j]))
      a <- as.vector(tcrossprod(x, y * kappa_val))
      can <- sum(exp(a + log(cpk))) / ngrid
      if (abs(can) < Inf) {
        den_mat[i, j] <- can
      }
    }
  }

  list(
    lat = lat_grid,
    long = long_grid,
    h = hw,
    den = den_mat
  )
}


#' Spherical density estimation
#'
#' @inheritParams sph_mean
# #' @param kamb logical. Whether to use the von Mises-Fisher kernel density estimation (`FALSE`) or Kamb's method (`TRUE`, the default).
#' @param FUN density estimation function if `kamb=TRUE`; one of [exponential_kamb()] (the default),
#'  [kamb_count], and [schmidt_count()].
#' @param n integer. Grid size. `128` by default.
#' @param sigma numeric. Radius for Kamb circle used for counting. 3 by default.
# #' @param vmf_hw numeric. Kernel bandwidth in degree.
# #' @param vmf_optimal character. Calculates an optimal kernel bandwidth
# #' using the cross-validation algorithm (`'cross'`) or the rule-of-thumb (`'rot'`)
# #' suggested by Garcia-Portugues (2013). Ignored when `vmf_hw` is specified.
#' @param weights (optional) numeric vector of length of `nrow(x)`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#' @inheritParams plot.Line
#' @param r numeric. radius of stereonet circle
#' @param ... arguments passed to [density_calc()]
#'
#' @name density
#' @aliases density.spherical density_spherical
#'
#' @examples
#' set.seed(20250411)
#' test <- rfb(100, mu = Line(120, 10), k = 5, A = diag(c(-1, 0, 1)))
#' density(x = test, n = 100, sigma = 3, weights = runif(100))
NULL

#' @rdname density
#' @exportS3Method stats::density
density.spherical <- function(x, ...) density_calc(x, ...)

#' @name density
#' @export
density_calc <- function(x,
                         # kamb = TRUE,
                         FUN = exponential_kamb,
                         n = 128L, sigma = 3,
                         # vmf_hw = NULL, vmf_optimal = c("cross", "rot"),
                         weights = NULL, upper.hem = FALSE, r = 1) {
  x_grid <- y_grid <- seq(-1, 1, length.out = n)

  grid <- expand.grid(x_grid, y_grid) |>
    as.matrix()
  dg <- density_grid(x, weights = weights, upper.hem = upper.hem, kamb = TRUE, FUN = FUN, sigma = sigma, ngrid = n, r = r)
  density_matrix <- matrix(dg$density, nrow = n, byrow = FALSE)
  dist_matrix <- grid[, 1]^2 + grid[, 2]^2

  # Create a logical mask where TRUE if outside the unit circle
  outside <- dist_matrix > r^2
  density_matrix[outside] <- NA
  res <- list(
    x = x_grid, y = y_grid,
    density = density_matrix
  )
  class(res) <- append(class(res), "sph_density")
  return(res)
}
