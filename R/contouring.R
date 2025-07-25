blank_grid_regular <- function(n, r = 1) {
  x_grid <- seq(-1, 1, length.out = n)
  y_grid <- seq(-1, 1, length.out = n)
  grid_schmidt <- expand.grid(x_grid, y_grid) |> as.matrix()

  # inside <- grid_schmidt[, 1]^2 + grid_schmidt[, 2]^2 <= r^
  # grid_schmidt <- grid_schmidt[inside, ]

  grid_cart <- .schmidt2cart(grid_schmidt[, 1], grid_schmidt[, 2])
  values <- rep(0, times = nrow(grid_cart))
  list(grid = grid_cart, density = values)
}

.schmdit2spherical <- function(x, y, r = 1) {
  tq <- sqrt(x^2 + y^2)

  # Recover inclination in radians
  inc_rad <- 2 * (pi / 4 - asin(tq / (r * sqrt(2))))

  # Recover azimuth in radians
  az_rad <- atan2(x, y) %% (2 * pi)

  # Return as matrix with named columns
  cbind(az_rad = az_rad, inc_rad = inc_rad)
}

.schmidt2cart <- function(x, y, r = 1) {
  res <- .schmdit2spherical(x, y, r) / DEG2RAD()
  lin2vec0(res[, 1], res[, 2])
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
#' @param azi,inc degree
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
  xyz_points <- lin2vec0(azi, inc)

  # cos_dist_matrix <- abs(grid_coords %*% t(xyz_points))

  # for(i in 1:ncol(cos_dist_matrix)){
  #   density_scale <- FUN(cos_dist_matrix[, i], sigma)
  # }

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

# count_points_fast <- function(azi, inc, FUN, sigma, n, weights, r) {
#   # Grid setup
#   grid <- blank_grid_regular(n = n, r = r)
#   grid_coords <- grid$grid
#
#   # Convert spherical to cartesian coordinates
#   xyz_points <- lin2vec0(azi, inc)
#
#   # Matrix multiplication: each row in grid_coords with all xyz_points
#   # t(grid_coords) %*% xyz_points computes all dot products efficiently
#   cos_dist_matrix <- abs(grid_coords %*% t(xyz_points)) # n_grid x n_points
#
#   # Apply FUN to entire matrix
#   density_scale <- FUN(cos_dist_matrix, sigma)
#
#   # Ensure weights are correctly recycled to n_points length
#   weights <- rep(weights, length.out = ncol(cos_dist_matrix))
#
#   # Compute weighted density sums across columns
#   weighted_density <- sweep(density_scale$count, 2, weights, "*")
#   density_sums <- rowSums(weighted_density)
#
#   # Compute final densities
#   densities <- (density_sums - 0.5) / density_scale$units
#   densities[densities < 0] <- 0
#
#   # Insert densities into grid
#   grid$density <- densities
#
#   grid
# }

#' Estimates point density
#'
#' Estimates point density of the given linear orientation measurements
#' Returns a regular (in lat-long space) grid of density estimates over a
#' hemispherical surface.
#'
#' @param x spherical object
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
  if (!is.spherical(x)) x <- to_spherical(x)

  if (!is.line(x)) x <- as.line(x)

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
      grid = Line(grid$azi, grid$inc) |> line2vec(),
      density = c(res$den)
    )
  }
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
    hw <- hw * DEG2RAD()
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



#' Calculate and plot densities in a stereonet
#'
#' Kamb counts and densities on the sphere
#'
#' @param x Object of class `"line"` or `"plane"` or `'spherical.density'` (for plotting only).
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param ngrid integer. Gridzise. 128 by default.
#' @param kamb logical. Whether to use the von Mises-Fisher kernel density estimation (`FALSE`) or Kamb's method (`TRUE`, the default).
#' @param FUN density estimation function if `kamb=TRUE`; one of [exponential_kamb()] (the default),
#'  [kamb_count], and [schmidt_count()].
#' @param sigma numeric. Radius for Kamb circle used for counting. 3 by default.
#' @param r numeric. radius of stereonet circle
#' @param weights (optional) numeric vector of length of `azi`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#' @param vmf_hw numeric. Kernel bandwidth in degree.
#' @param vmf_optimal character. Calculates an optimal kernel bandwidth
#' using the cross-validation algorithm (`'cross'`) or the rule-of-thumb (`'rot'`)
#' suggested by Garcia-Portugues (2013). Ignored when `vmf_hw` is specified.
#' @param nlevels integer. Number of contour levels for plotting
#' @param type character. Type of plot: `'contour'` for contour lines,
#' `'contour_filled'` for filled contours, or `'image'` for a raster image.
#' @param add logical. Whether the contours should be added to an existing plot.
#' @param col colour(s) for the contour lines drawn. If `NULL`, lines are color based on `col.palette`.
#' @param col.palette a color palette function to be used to assign colors in the plot.
#' @param col.params list. Arguments passed to `col.palette`
#' @param ... optional parameters passed to [image()] or [contour()].
#'
#' @name stereo_density
#'
#' @seealso [count_points()]
#'
#' @returns list containing the stereographic x and coordinates of of the grid,
#' the counts, and the density.
#'
#' @examples
#' set.seed(20250411)
#' test <- rfb(100, mu = Line(120, 10), k = 5, A = diag(c(-1, 0, 1)))
#' test_densities <- spherical_density(x = test, ngrid = 100, sigma = 3, weights = runif(100))
#'
#' stereo_density(test_densities, type = "image", add = FALSE)
#' stereo_point(test, col = "lightgrey", pch = 21)
#'
#' stereo_density(test,
#'   type = "contour_filled", add = FALSE,
#'   col.params = list(direction = -1, begin = .05, end = .95, alpha = .75)
#' )
#' stereo_point(test, col = "black", pch = 21)
#'
#' # complete example:
#' par(mfrow = c(2, 1))
#' wp <- 6 / ifelse(is.na(example_planes$quality), 6, example_planes$quality)
#' my_planes <- Plane(example_planes$dipdir, example_planes$dip)
#' fabric_p <- or_shape_params(my_planes)$Vollmer["D"]
#' my_planes_eig <- or_eigen(my_planes)
#'
#' stereoplot(guides = TRUE, col = "grey96")
#' stereo_point(my_planes, col = "grey", pch = 16, cex = .5)
#' stereo_density(my_planes, type = "contour", add = TRUE, weights = wp)
#' stereo_point(as.plane(my_planes_eig$vectors[3, ]), col = "black", pch = 16)
#' stereo_greatcircle(as.plane(my_planes_eig$vectors[3, ]), col = "black", pch = 16)
#' title(
#'   main = "Planes",
#'   sub = paste0(
#'     "N: ", nrow(my_planes), " | Fabric strength: ", round(fabric_p, 2),
#'     "\nLambert equal area, lower hemisphere projection"
#'   )
#' )
#'
#' my_lines <- Line(example_lines$trend, example_lines$plunge)
#' wl <- 6 / ifelse(is.na(example_lines$quality), 6, example_lines$quality)
#' fabric_l <- or_shape_params(my_lines)$Vollmer["D"]
#'
#' stereoplot(guides = TRUE, col = "grey96")
#' stereo_point(my_lines, col = "grey", pch = 16, cex = .5)
#' stereo_density(my_lines, type = "contour", add = TRUE, weights = wl)
#' stereo_point(v_mean(my_lines, w = wl), col = "black", pch = 16)
#' title(
#'   main = "Lines",
#'   sub = paste0(
#'     "N: ", nrow(my_lines), " | Fabric strength: ", round(fabric_l, 2),
#'     "\nLambert equal area, lower hemisphere projection"
#'   )
#' )
NULL

#' @rdname stereo_density
#' @export
spherical_density <- function(x,
                              kamb = TRUE, FUN = exponential_kamb,
                              ngrid = 128L, sigma = 3,
                              vmf_hw = NULL, vmf_optimal = c("cross", "rot"),
                              weights = NULL, upper.hem = FALSE, r = 1) {
  x_grid <- y_grid <- seq(-1, 1, length.out = ngrid)

  grid <- expand.grid(x_grid, y_grid) |> as.matrix()

  dg <- if (kamb) {
    density_grid(x, weights = weights, upper.hem = upper.hem, kamb = TRUE, FUN = FUN, sigma = sigma, ngrid = ngrid, r = r)
  } else {
    stop("vmf not supported at the moment")
    # density_grid(x, weights = NULL, upper.hem = upper.hem, kamb = FALSE, ngrid = ngrid, hw = vmf_hw, kernel_method = vmf_optimal)
  }


  # if(!kamb){
  #   grd_lines <- dg$grid |> vec2line()
  #   grid <- .schmidt_crds(deg2rad(grd_lines[, 1]), deg2rad(grd_lines[, 2]), r = r)
  #
  #   grid <- cbind(grid, dg$density)
  #   grid <- grid[order(grid[, 1], grid[, 2]), ]
  #
  #   x_grid <- unique(grid[, 1])
  #   y_grid <- unique(grid[, 2])
  #   density_matrix <- matrix(grid[, 3], nrow = ngrid, byrow = FALSE)
  #
  # } else{
  density_matrix <- matrix(dg$density, nrow = ngrid, byrow = FALSE)
  # }

  dist_matrix <- grid[, 1]^2 + grid[, 2]^2


  # Create a logical mask where TRUE if outside the unit circle
  outside <- dist_matrix > r^2

  density_matrix[outside] <- NA

  res <- list(
    x = x_grid, y = y_grid,
    density = density_matrix
  )
  class(res) <- append(class(res), "spherical.density")
  return(res)
}

# #' @rdname stereo_density
# #' @export
# projected_density <- function(x, ngrid = 128L, sigma = 3, weights = NULL, upper.hem = FALSE, r = 1) {
#   if (is.plane(x) | is.fault(x)) {
#     x[, 1] <- 180 + x[, 1]
#     x[, 2] <- 90 - x[, 2]
#   }
#
#   crds <- stereo_coords(
#     x[, 1],
#     x[, 2],
#     upper.hem,
#     earea = TRUE
#   )
#
#
#   # prepare the grid
#   x_grid <- y_grid <- seq(-1, 1, length.out = ngrid)
#   grid <- expand.grid(x_grid, y_grid) |> as.matrix()
#
#   # Compute distance from origin
#   dist_matrix <- sqrt(grid[, 1]^2 + grid[, 2]^2)
#
#   # Create a logical mask where TRUE if outside the unit circle
#   mask_outside <- dist_matrix > r
#
#   N <- nrow(crds)
#
#   if (is.null(weights)) {
#     weights <- rep(1, N)
#   }
#
#
#   # f <- sigma^2 / (sigma^2 + N)
#   counter_radius <- sigma * r / sqrt(N + sigma^2)
#
#   counts <- .count_points_within_radius(crds, grid, counter_radius)
#
#
#   counts_matrix <- matrix(counts, nrow = n, byrow = FALSE)
#   counts_matrix[mask_outside] <- NA # replace cells outside of unit circle with NA
#
#   density_matrix <- counts_matrix / max(counts_matrix, na.rm = TRUE)
#
#   res <- list(
#     x = x_grid, y = y_grid,
#     grid = grid,
#     counts = counts_matrix,
#     density = density_matrix
#   )
#   class(res) <- append(class(res), "spherical.density")
#   return(res)
# }

#' @rdname stereo_density
#' @export
stereo_density <- function(x, kamb = TRUE, FUN = exponential_kamb, ngrid = 128L, sigma = 3,
                           vmf_hw = NULL, vmf_optimal = c("cross", "rot"),
                           weights = NULL, upper.hem = FALSE, r = 1,
                           type = c("contour", "contour_filled", "image"), nlevels = 10L,
                           col.palette = viridis, col = NULL, add = TRUE, col.params = list(),
                           ...) {
  type <- match.arg(type)

  if (inherits(x, "spherical.density")) {
    d <- x
  } else {
    d <- spherical_density(x,
      ngrid = ngrid,
      kamb = kamb,
      upper.hem = upper.hem, r = r,
      sigma = sigma, FUN = FUN, weights = weights,
      vmf_hw = vmf_hw, vmf_optimal = vmf_optimal
    )
  }

  densities <- d$density

  if (type == "image") {
    if (!add) {
      stereoplot(guides = FALSE)
      add <- TRUE
    }

    col.params <- append(list(n = nlevels), col.params)
    col <- do.call(col.palette, col.params)

    graphics::image(
      x = d$x, y = d$y,
      z = densities,
      col = col,
      asp = 1,
      axes = FALSE,
      frame.plot = FALSE,
      add = add,
      ...
    )
  } else if (type == "contour") {
    if (!add) {
      stereoplot(guides = FALSE)
      add <- TRUE
    }

    if (is.null(col)) {
      levels <- pretty(range(densities, na.rm = TRUE), nlevels)
      col.params <- append(list(n = length(levels) - 1), col.params)
      col <- do.call(col.palette, col.params)
    }

    graphics::contour(
      x = d$x, y = d$y,
      z = densities,
      levels = levels,
      col = col,
      asp = 1,
      axes = FALSE,
      frame.plot = FALSE,
      add = add,
      ...
    )
  } else {
    if (!add) stereoplot(guides = FALSE)

    levels <- pretty(range(densities, na.rm = TRUE), nlevels)
    col.params <- append(list(n = length(levels) - 1), col.params)
    col <- do.call(col.palette, col.params)

    graphics::.filled.contour(
      x = d$x, y = d$y,
      z = densities,
      levels = levels,
      col = col
    )
  }
}


# .count_points_within_radius <- function(A, B, r) {
#   # Compute squared distances efficiently
#   d2 <- (outer(A[, 1], B[, 1], "-"))^2 +
#     (outer(A[, 2], B[, 2], "-"))^2
#
#   # Logical matrix where TRUE if within radius
#   within <- d2 <= r^2
#
#   # Count per column (i.e., per grid point)
#   counts <- colSums(within)
#
#   return(counts)
# }

#
# .fix_symm <- function(x) {
#   x[x[, 3] < 0, ] <- v_antipode(x[x[, 3] < 0, ])
#   x
# }
