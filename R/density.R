# # Simulation of random values from rotationally symmetric distributions
# x <- rvmf(n = 200, mu = Line(120, 30), k = 15) |> to_vec()
#
# # MLE of (hyper-)spherical rotational symmetric distributions
# Directional::kent.mle(x)
# Directional::vmf.mle(x)
#
# # Tuning of the bandwidth parameter in the von Mises-Fisher kernel
# h <- Directional::vmfkde.tune(x)[1] # estimate
#
# # von Mises-Fisher kernel density estimation
# Directional::vmf.kde(x, h = h, thumb = "none")
# Directional::vmf.kde(x, h = h, thumb = "rot")
#
# # von Mises-Fisher kernel density estimate
# x_den <- Directional::vmf.kerncontour(Directional::euclid.inv(x), thumb = "none", full = TRUE, ngrid = 100, den.ret = TRUE)
#
# # Contour plot of spherical data using a von Mises-Fisher kernel density estimate
# filled.contour(
#   x_den$lat, x_den$long, x_den$den,
#   nlevels = 200, color.palette = function(n) scico::scico(n),
#   plot.axes = {
#     axis(1, col = "black", cex.axis = 1.2)
#     axis(2, col = "black", cex.axis = 1.2)
#     contour(x_den$lat, x_den$long, x_den$den,
#             col = "black", nlevels = 10,
#             labcex = 0.8, lwd = 1.5, add = TRUE)
#     }, key.axes = {
#       axis(4, col = "black", cex.axis = 1.2)
#       },
#   xlab = "Azimuth (=Latitude)", ylab = "Plunge (=Longitude)", cex.lab = 1.2)
#
#
# vmf_kde <- function(x, ngrid) {
#   class = class(x)
#   x <- to_vec(x)
#   n <- nrow(x)
#
#   h <- as.numeric(Directional::vmfkde.tune(x, low = 0.1, up = 1)[1])
#
#   # x1 <- seq(0, 180, length = ngrid) # lat
#   # x2 <- seq(0, 360, length = grid) # long
#
#   x1 <- seq(-90, 90, length = ngrid) # lat
#   x2 <- seq(-180, 180, length = ngrid) # long
#
#   cpk <- 1 / ((h^2)^0.5 * (2 * pi)^1.5 * besselI(1 / h^2, 0.5))
#   mat <- matrix(nrow = ngrid, ncol = ngrid)
#   for (i in 1:ngrid) {
#     for (j in 1:ngrid) {
#       y <- geo2vec(c(x1[i], x2[j]))
#       a <- as.vector(tcrossprod(x, y / h^2))
#       can <- sum(exp(a + log(cpk))) / ngrid
#       if (abs(can) < Inf) {
#         mat[i, j] <- can
#       }
#     }
#   }
#
#   # geographical to spherical coordinates
#   sph <- geo2vec(cbind(x1-90, x2-180)) |> to_spherical(class)
#
#   return(list(azi = sph[, 1], inc = sph[, 2], h = h, den = mat))
# }
#
# vmf_kde_grid <- function(x, ngrid = 100, upper.hem = FALSE) {
#   # Translate to (0,180) and (0,360)
#   # x[, 1] <- x[, 1]  # lat
#   # x[, 2] <- x[, 2]  # long
#   # x <- to_vec(x) |> Directional::euclid.inv()
#
#   res <- vmf_kde(x, ngrid = ngrid)
#
#   sph <- cbind(res$azi, res$inc) |> to_vec() |> vec2line()
#
#   # spherical to stereographic
#   crds <- stereo_coords(sph[, 1], sph[, 2], upper.hem)
#
#
#   grid <- expand.grid(x = crds[, 1], y = crds[, 2])
#   grid$density <- c(res$den)
#   grid
# }
# rvmf(n = 200, mu = Line(120, 30), k = 15)  |>
# vmf_kde_grid() |>
# ggplot(aes(x, y, color = density)) +
# geom_point() +
# scico::scale_color_scico(palette = "bilbao")


#' Density distribution of vectors
#'
#' Calculates density distribution of vectors using the modified Kamb contouring with exponential smoothing.
#'
#' @param x object of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param sigma numeric. If `NULL` sigma is calculated automatically. Default `NULL`
#' @param sigmanorm logical. If `TRUE` counting is normalized to sigma multiples. Default `TRUE`
#' @param trimzero logical. If `TRUE` zero contour is not drawn. Default `FALSE`
#' @returns list
#' @examples
#' set.seed(20250411)
#' x <- rvmf(n = 200, mu = Line(120, 30), k = 10)
#' calculate_density(x)
# calculate_density <- function(x, sigma = NULL, sigma.norm = TRUE, trimzero = TRUE, ngrid = 3000, grid.type = c("sfs", "gss", "rot")) {
#   x <- Vec3(x) |> unclass()
#   # parse options
#   grid.type <- match.arg(grid.type)
#   n <- nrow(x)
#
#   if (is.null(sigma)) {
#     # k = estimate_k(x)
#     # sigma = sqrt(2 * n / (k - 2))
#     # Totally empirical as estimate_k is problematic
#     sigma <- sqrt(2 * n / (log(n) - 2)) / 3
#     k <- 2 * (1 + n / sigma^2)
#   } else {
#     k <- 2 * (1 + n / sigma^2)
#   }
#   # method = kwargs.get("method", "exp_kamb")
#
#   # do calc
#   scale <- sqrt(n * (k / 2 - 1) / k^2)
#
#   grid <- runif.spherical(n = ngrid, method = grid.type) |> unclass()
#
#   cnt <- exp(k * abs(grid %*% t(x)) - 1)
#   x2 <- colSums(cnt) / scale
#
#   if (sigma.norm) {
#     x2 <- x2 / sigma
#   }
#
#   x2[which(x2 < 0)] <- 0
#   # x2 <- abs(x)
#
#   if (trimzero) {
#     x2[which(x2 == 0)] <- .Machine$double.xmin
#   }
#
#   # density_params
#   list(
#     grid = grid,
#     density = x2,
#     k = k,
#     scale = scale,
#     sigma = sigma,
#     sigma.norm = sigma.norm
#   )
# }

#' Contour lines of density in stereographic projection
#'
#' plots contour lines of density in stereographic projection.
#'
#' @inheritParams calculate_density
#' @param ... arguments passed to [graphics::contour()]
#' @importFrom graphics filled.contour
#' @importFrom stats xtabs
# stereo_density_contour <- function(x, sigma = NULL, sigma.norm = TRUE, trimzero = TRUE, ngrid = 100, grid.type = c("gss", "sfs"), ...) {
#   grid.type <- match.arg(grid.type)
#
#   densgrd <- calculate_density(
#     x,
#     sigma = sigma, sigma.norm = sigma.norm, trimzero = trimzero, ngrid = ngrid, grid.type = grid.type
#   )
#
#   XY <- project_data(
#     x,
#     x = densgrd$grid[, 1],
#     y = densgrd$grid[, 2],
#     z = densgrd$grid[, 3],
#     clip.inside = FALSE, upper.hem = TRUE, rotate.data = FALSE
#   ) / (3 / 2)
#
#   d_stereo <- cbind(XY, densgrd$density)
#   colnames(d_stereo) <- c("x", "y", "d")
#
#   m <- xtabs(d ~ x + y, data = d_stereo)
#
#   filled.contour(as.numeric(rownames(m)), as.numeric(colnames(m)), m, nlevels = 20, color.palette = viridis::magma)
# }

# blank_grid <- function(n = 3000, ...) {
#   grid <- v_unif(NULL, n = n, ...)
#   values <- rep(0, times = n)
#   list(grid = grid, density = values)
# }


#' @importFrom dplyr near
# project_data0 <- function(self, x, y, z) {
#   # Equal-area projection
#   d <- sqrt(x * x + y * y + z * z)
#   if (any(d == 0)) {
#     return(cbind(NA, NA))
#   } else {
#     x <- x / d
#     y <- y / d
#     z <- z / d
#
#     # z[np.isclose(1 + z, np.zeros_like(z))] = 1e-6 - 1
#     z[which(near(1 + z, rep(0, length(z))))] <- 1e6 - 1
#
#     sqz <- sqrt(1 / (1 + z))
#     return(cbind(y * sqz, x * sqz))
#   }
# }

# project_data <- function(self, x, y, z, clip.inside = TRUE, upper.hem = FALSE, rotate.data = FALSE) {
#   if (rotate.data) {
#     xyz <- vresultant(x) %>% t(cbind(x, y, z))
#     x <- xyz[, 1]
#     y <- xyz[, 2]
#     z <- xyz[, 3]
#   }
#   if (upper.hem) {
#     XY <- project_data0(self, -x, -y, -z)
#     X <- XY[, 1]
#     Y <- XY[, 2]
#     if (clip.inside) {
#       outside <- X * X + Y * Y > 1
#       X[which(outside)] <- NA
#       Y[which(outside)] <- NA
#     }
#     return(cbind(-X, -Y))
#   } else {
#     XY <- project_data0(self, x, y, z)
#     X <- XY[, 1]
#     Y <- XY[, 2]
#
#     if (clip_inside) {
#       outside <- X * X + Y * Y > 1
#       X[which(outside)] <- NA
#       Y[which(outside)] <- NA
#     }
#     return(cbind(x = X, y = Y))
#   }
# }



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
  res <- .schmdit2spherical(x, y, r)  |> rad2deg()
  Line(res[, 1], res[, 2]) |> Vec3() |> unclass()
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
  xyz_points <- Line(azi, inc) |> Vec3() |> unclass()

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
  # if (!is.spherical(x)) x <- to_spherical(x)

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
#' @param x object of class `"Vec3()"`, `"Line()"`, or `"Plane()"`.
#' @param kamb logical. Whether to use the von Mises-Fisher kernel density estimation (`FALSE`) or Kamb's method (`TRUE`, the default).
#' @param FUN density estimation function if `kamb=TRUE`; one of [exponential_kamb()] (the default),
#'  [kamb_count], and [schmidt_count()].
#' @param ngrid integer. Gridzise. 128 by default.
#' @param sigma numeric. Radius for Kamb circle used for counting. 3 by default.
#' @param vmf_hw numeric. Kernel bandwidth in degree.
#' @param vmf_optimal character. Calculates an optimal kernel bandwidth
#' using the cross-validation algorithm (`'cross'`) or the rule-of-thumb (`'rot'`)
#' suggested by Garcia-Portugues (2013). Ignored when `vmf_hw` is specified.
#' @param weights (optional) numeric vector of length of `azi`.
#' The relative weight to be applied to each input measurement. The array
#' will be normalized to sum to 1, so absolute value of the `weights` do not
#' affect the result. Defaults to `NULL`
#' @inheritParams plot.spherical
#' @param r numeric. radius of stereonet circle
#'
#' @rdname stereo_density
#' @export
#'
#' @examples
#' set.seed(20250411)
#' test <- rfb(100, mu = Line(120, 10), k = 5, A = diag(c(-1, 0, 1)))
#' density(x = test, ngrid = 100, sigma = 3, weights = runif(100))
density.spherical <- function(x,
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

