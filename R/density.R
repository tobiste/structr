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


#'
#'
#' @title Density distribution of vectors
#'
#' Calculates density distribution of vectors using the modified Kamb contouring with exponential smoothing.
#'
#' @param x description
#' @param sigma numeric. If `NULL` sigma is calculated automatically. Default `NULL`
#' @param sigmanorm logical. If `TRUE` counting is normalized to sigma multiples. Default `TRUE`
#' @param trimzero logical. If `TRUE` zero contour is not drawn. Default `FALSE`
#' @returns description
#' @examples
#' @examples
#' x <- rvmf(n = 200, mu = Line(120, 30), k = 10)
#' calculate_density(x)
calculate_density <- function(x, sigma = NULL, sigma.norm = TRUE, trimzero = TRUE, ngrid = 3000, grid.type = c("sfs", "gss", "rot")) {
  x <- to_vec(x)
  # parse options
  grid.type <- match.arg(grid.type)
  n <- nrow(x)

  if (is.null(sigma)) {
    # k = estimate_k(x)
    # sigma = sqrt(2 * n / (k - 2))
    # Totally empirical as estimate_k is problematic
    sigma <- sqrt(2 * n / (log(n) - 2)) / 3
    k <- 2 * (1.0 + n / sigma^2)
  } else {
    k <- 2 * (1.0 + n / sigma^2)
  }
  # method = kwargs.get("method", "exp_kamb")

  # do calc
  scale <- sqrt(n * (k / 2 - 1) / k^2)

  grid <- v_unif(NULL, ngrid, method = grid.type)

  cnt <- exp(k * abs(grid %*% t(x)) - 1)
  x2 <- colSums(cnt) / scale

  if (sigma.norm) {
    x2 <- x2 / sigma
  }

  x2[which(x2 < 0)] <- 0
  # x2 <- abs(x)

  if (trimzero) {
    x2[which(x2 == 0)] <- .Machine$double.xmin
  }

  # density_params
  list(
    grid = grid,
    density = x2,
    k = k,
    scale = scale,
    sigma = sigma,
    sigma.norm = sigma.norm
  )
}

#' Contour lines of density in stereographic projection
#'
#' plots contour lines of density in stereographic projection.
#'
#' @inheritParams calculate_density
#' @param ... arguments passed to [graphics::contour()]
#' @examples
#' x <- rvmf(n = 200, mu = Line(120, 30), k = 10)
#' stereoplot()
#' stereo_density_contour(x)
#' stereo_point(x)
stereo_density_contour <- function(x, sigma = NULL, sigma.norm = TRUE, trimzero = TRUE, ngrid = 100, grid.type = c("gss", "sfs"), ...) {
  # stereoplot()
  densgrd <- calculate_density(
    x,
    sigma = sigma, sigma.norm = sigma.norm, trimzero = trimzero, ngrid = ngrid, grid.type = grid.type
  )

  XY <- project_data(
    x,
    x = densgrd$grid[, 1],
    y = densgrd$grid[, 2],
    z = densgrd$grid[, 3],
    clip.inside = FALSE, upper.hem = TRUE, rotate.data = FALSE
  ) / (3 / 2)

  d_stereo <- cbind(XY, densgrd$density)
  colnames(d_stereo) <- c("x", "y", "d")

  m <- xtabs(d ~ x + y, data = d_stereo)

  filled.contour(as.numeric(rownames(m)), as.numeric(colnames(m)), m, nlevels = 20, color.palette = viridis::magma)
}
#'
#'
blank_grid <- function(n = 3000, ...) {
  grid <- v_unif(NULL, n = n, ...)
  values <- rep(0, times = n)
  list(grid = grid, density = values)
}
#'
#'
project_data0 <- function(self, x, y, z) {
  # Equal-area projection
  d <- sqrt(x * x + y * y + z * z)
  if (any(d == 0)) {
    return(cbind(NA, NA))
  } else {
    x <- x / d
    y <- y / d
    z <- z / d

    # z[np.isclose(1 + z, np.zeros_like(z))] = 1e-6 - 1
    z[which(near(1 + z, rep(0, length(z))))] <- 1e6 - 1

    sqz <- sqrt(1 / (1 + z))
    return(cbind(y * sqz, x * sqz))
  }
}
#'
#'
#'
project_data <- function(self, x, y, z, clip.inside = TRUE, upper.hem = FALSE, rotate.data = FALSE) {
  if (rotate.data) {
    xyz <- vresultant(x) %>% t(cbind(x, y, z))
    x <- xyz[, 1]
    y <- xyz[, 2]
    z <- xyz[, 3]
  }
  if (upper.hem) {
    XY <- project_data0(self, -x, -y, -z)
    X <- XY[, 1]
    Y <- XY[, 2]
    if (clip.inside) {
      outside <- X * X + Y * Y > 1
      X[which(outside)] <- NA
      Y[which(outside)] <- NA
    }
    return(cbind(-X, -Y))
  } else {
    XY <- project_data0(self, x, y, z)
    X <- XY[, 1]
    Y <- XY[, 2]

    if (clip_inside) {
      outside <- X * X + Y * Y > 1
      X[which(outside)] <- NA
      Y[which(outside)] <- NA
    }
    return(cbind(x = X, y = Y))
  }
}
