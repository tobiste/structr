#' \deqn{R_f/\phi}
#'
#' @param Rs numeric.
#' @param Ri numeric. Initial or pre-deformation aspect ratio of an elliptical object
#' @param theta numeric. Fluctuation angle of undeformed marker w.r.t. computed principal direction of strain (in degrees)
#'
#' @returns list. `phi`: the final angle (in degrees) between long axis and principal strain direction, `Rf`: final axial ratio
#' @export
#'
#' @examples
#' #' data(ramsay)
#' Rphi(ramsay[, "R"], theta = ramsay[, "phi"])
Rphi <- function(Rs, Ri = 1, theta) {
  theta <- theta * pi / 180
  a <- 2 * Rs * (Ri^2 - 1) * sin(2 * theta)
  b <- (Ri^2 + 1) * (Rs - 1) + (Ri^2 - 1) * (Rs^2 + 1) * cos(2 * theta)
  tan_2phi <- a / b
  phi <- atan(tan_2phi) / 2

  c <- tan(phi)^2 * (1 + Ri^2 * tan(theta)^2) - Rs^2 * (tan(theta)^2 + Ri^2)
  d <- Rs^2 * tan(phi)^2 * (tan(theta)^2 + Ri^2) - (1 + Ri^2 * tan(theta)^2)
  Rf <- sqrt(c / d)

  return(list("phi" = phi * pi / 180, "Rf" = Rf))
}


#' Mean strain ellipse
#'
#' determines the shape and orientation of the strain ellipse by using
#' deformed elliptical objects as strain markers. The algorithm is based on the
#' mean shape matrix and its eigenvalues (Shimamoto and Ikeda, 1976).
#'
#' @param r numeric. Aspect ratio of deformed object (long axis / short axis)
#' @param phi numeric. Orientation of long axis of deformed object (in degrees)
#' @param boot logical. Whether a 95% confidence interval from on bootstrapping should be calculated. `TRUE` by default.
#' @param iter integer. Number of bootstrap resamples (`1000` by default). Ignored when `boot = FALSE`.
#' @param boot.values logical. Whether the bootstrapped R and phi values should be added to the output. `FALSE` by default.
#'
#' @returns list. `R` gives the mean aspect ratio of the strain ellipse,
#' and `phi` gives the orientation of its long axis. If `boot=TRUE`, then
#' the bootstrapped 95% confidence interval for the mean aspect ratio (`R_CI`)
#' and for its orientation (`phi_CI`) are added. If `boot.values=TRUE`, a matrix
#' containing the bootstrapped R and phi values are added.
#'
#' @export
#'
#' @references Shimamoto, T., Ikeda, Y., 1976. A simple algebraic method for
#' strain estimation from ellipsoidal objects. Tectonophysics 36, 315–337.
#' doi: 10.1016/0040-1951(76)90107-4
#'
#' @examples
#' set.seed(20250411)
#' data(ramsay)
#' mean_strain_ellipse(ramsay[, "R"], ramsay[, "phi"])
mean_strain_ellipse <- function(r, phi, boot = TRUE, resamples = 1000, boot.values = FALSE) {
  res <- mean_strain_ellipse0(r, phi)
  if (boot) {
    rphi_boot <- vapply(1:resamples, function(i) {
      rsmpl <- sample(seq_along(r), replace = TRUE)
      mse <- mean_strain_ellipse0(r[rsmpl], phi[rsmpl])
      unlist(mse)
    }, FUN.VALUE = numeric(2)) |>
      t()

    res$R_CI <- quantile(rphi_boot[, 1], probs = c(.025, .975)) |> unname()
    res$phi_CI <- quantile(rphi_boot[, 2], probs = c(.025, .975)) |> unname()

    if (boot.values) res$boot <- rphi_boot
  }
  res
}

mean_strain_ellipse0 <- function(r, phi) {
  n <- length(r)
  stopifnot(n == length(phi))
  phi <- phi * pi / 180
  cosp2 <- cos(phi)^2
  sinp2 <- sin(phi)^2
  f <- cosp2 / r + r * sinp2
  g <- sinp2 / r + r * cosp2
  h <- ((1 / r) - r) * sin(phi) * cos(phi)

  if (n > 1) {
    f <- mean(f)
    g <- mean(g)
    h <- mean(h)
  }

  sm <- matrix(c(f, h, h, g), nrow = 2, ncol = 2)


  theta <- (atan((2 * h / (f - g))) / 2) * 180 / pi

  sm_eig <- eigen(sm, symmetric = TRUE, only.values = TRUE)
  pa <- sqrt(sm_eig$values)

  list(
    R = pa[1] / pa[2],
    phi = theta
  )
}




#' Calculates densities for fabric and strain data
#'
#' Densities in hyperbaloidal projections of geological fabric and
#' finite strain data density calculations done on the unit hyperbaloid
#' (Vollmer, 2018). Options are given for equidistant (Elliott), equal-area,
#' stereographic, orthographic, exponential, and radial projections, as polar
#' azimuthal or cylindrical (cartesian, RfPhi-type) plots.
#'
#' @inheritParams mean_strain_ellipse
#' @param rmax maximum R value (if `NULL`, computed automatically)
#' @param kappa smoothing parameter
#' @param nnodes grid resolution. higher is more accurate but slower, 30 is good, but 50 is recommended for final plots, default = 50.
# #' @param angfmt character. Angle format of `phi`: `'deg'` for degree (the default), `'rad'` for radians, and `'grd'` for gradians.
#' @param proj character. Projection,  `'eqd'` for equidistant (Elliot plot),
#' `'eqa'` for equal-area, `'stg'` for stereographic, `'ort'` for orthographic,
#' `'gno'` for gnomonic, `'lin'` for exponential (linear R), `'rdl'` for radial,
#' `'rfp'` for Rf/phi (cylindrical instead of polar).
#' @param normalize logical.
#'
#' @return list with:
#'   - x, y: vectors of the projection coordinates for the density grid
#'   - z: density matrix
#'   - points: projected data points
#'   - frame: plot frame (circle or square)
#'
#' @details
#' The data must be in a comma delimited csv text file with one (R, phi) pair
#' per line, where R = strain ratio (max/min), phi = orientation of long (max)
#' axis from x. Contours are equally spaced over the probability density
#' distribution.
#'
#' @references
#' Vollmer, F.W., 2018. Automatic contouring of geological fabric and finite
#' strain data on the unit hyperboloid. Computers & Geosciences,
#' https://doi.org/10.1016/j.cageo.2018.03.006
#'
#' @export
#'
#' @seealso [Rphi_plot()], [Rphi_polar_plot()]
#'
#' @examples
#' data(ramsay)
#' out <- hypercontour(ramsay[, "R"], ramsay[, "phi"], proj = "rfp")
#' head(out$z)
hypercontour <- function(r, phi,
                         # angfmt = c("deg", "rad", "grd"),
                         proj = c("eqd", "eqa", "stg", "ort", "gno", "lin", "rdl", "rfp"),
                         normalize = TRUE,
                         rmax = NULL, kappa = 40, nnodes = 50L) {
  # angfmt <- match.arg(angfmt)
  proj <- match.arg(proj)


  frame <- .drawFrame(TRUE, proj)

  if (is.null(rmax)) {
    rmax <- ceiling(max(r)) + 1
  }

  # convert angles
  # f <- switch(angfmt,
  #   "deg" = pi / 180, # degrees
  #   "grd" = pi / 200, # gradians
  #   "rad" = pi
  # ) # radians
  f <- pi / 180
  phi <- phi * f

  points <- .Rphi2xy(r, phi, rmax = rmax, proj = proj)
  grd <- gridHyper(cbind(R = r, phi = phi), rmax = rmax, kappa = kappa, nnodes = nnodes, normalize = normalize, proj = proj)

  x <- seq(-1, 1, length.out = nnodes)
  y <- seq(-1, 1, length.out = nnodes)
  z <- grd

  if (proj == "rfp") {
    # xy <- .xy2Rphi(x, y, rmax=rmax, proj = 'rfp')
    x <- (x * pi / 2) / f
    y <- .scale(y, to = c(0, rmax))
    #
    frame[, "x"] <- (frame[, "x"] * pi / 2) / f
    frame[, "y"] <- .scale(frame[, "y"], to = c(0, rmax))
    #
    # points <- .xy2Rphi(points[, 1], points[, 2], rmax, proj='rfp')

    points <- cbind(x = phi / f, y = r)
    # points[, 'phi'] <- points[, 'phi'] / f
  } else {
    xy <- expand.grid(x, y)
    r <- 1
    r2 <- r * r
    mask <- xy[, 1]^2 + xy[, 2]^2 > r2
    z[mask] <- NA
  }

  return(list(
    x = x,
    y = y,
    z = z,
    frame = frame,
    points = points,
    proj = proj,
    rmax = rmax
  ))
}

.R2zeta <- function(r, proj) {
  switch(as.character(proj),
    "eqd" = log(r),
    "eqa" = sqrt(r) - 1 / sqrt(r),
    "stg" = {
      t <- sqrt(r)
      s <- 1 / t
      2 * (t - s) / (t + s)
    },
    "ort" = 0.5 * (r - 1 / r),
    "gno" = {
      t <- r * r
      (t - 1) / (t + 1)
    },
    "lin" = r - 1,
    "rdl" = 0.5 * (r + 1 / r) - 1,
    "rfp" = r,
    stop("Unknown proj in rToZeta")
  )
}


.zeta2R <- function(z, proj) {
  switch(as.character(proj),
    "eqd" = exp(z),
    "eqa" = {
      t <- z + sqrt(z^2 + 4)
      (t^2) / 4
    },
    "stg" = {
      t <- z / 2
      (1 + t) / (1 - t)
    },
    "ort" = z + sqrt(z^2 + 1),
    "gno" = {
      t <- 0
      if (z < 0.99) t <- sqrt((1 + z) / (1 - z))
      if (t < 50.001) t else 0
    },
    "lin" = z + 1,
    "rdl" = {
      t <- z + 1
      t + sqrt(t^2 - 1)
    },
    "rfp" = z,
    stop("Unknown proj in zetaToR")
  )
}

# rPhiToXY - projects R, phi to cartesian coordinates of unit hyperbaloidal
# projection. Maps to [-1..-1, +1..+1] to overlie unit image.
.Rphi2xy <- function(r, phi, rmax, proj) {
  # stopifnot(is.logical(rfp))
  z <- .R2zeta(r, proj)
  zm <- .R2zeta(rmax, proj)
  s <- z / zm
  if (proj == "rfp") {
    p <- phi
    p <- ifelse(p < -pi / 2, p + pi, p)
    p <- ifelse(p >= pi / 2, p - pi, p)
    x <- 2 * p / pi
    y <- 2 * s - 1
  } else {
    x <- s * cos(2 * phi)
    y <- s * sin(2 * phi)
  }

  cbind(x = x, y = y)
}


# xYToRPhi - back projects cartesian coordinates of hyperbaloidal projection.
# Not scaled from unit plot.
.xy2Rphi <- function(x, y, rmax, proj) {
  zm <- .R2zeta(rmax, proj)
  if (proj == "rfp") {
    z <- (y + zm) / 2
    r <- .zeta2R(z, proj)
    phi <- x * (0.5 * pi / zm)
    if (phi < -pi / 2) {
      phi <- phi + pi
    } else if (phi >= pi / 2) phi <- phi - pi
  } else { # polar
    t <- sqrt(x * x + y * y)
    r <- .zeta2R(t, proj)
    phi <- atan2(y, x) / 2
  }
  cbind(R = r, phi = phi)
}

# rhoPsiToH - set as a hyperbolic position vector from rho and psi. For strain
# ellipses: rho = ln(R), psi = 2 phi. Ref: Yamaji, 2008. }
.rhopsi2H <- function(rho, psi) {
  s <- sinh(rho)
  x <- cosh(rho)
  y <- s * cos(psi)
  z <- s * sin(psi)
  cbind(x = x, y = y, z = z)
}

#  converts R, phi to hyperbaloidal point.
.Rphi2H <- function(r, phi) {
  rho <- ifelse(r < 1, 0, log(r))

  psi <- 2 * phi
  .rhopsi2H(rho, psi)
}

# dotH - hyperbolic inner product. Ref: Yamaji, 2008, eqn 4.
.dotH <- function(hn, H) {
  (-hn[, 1] * H[, 1]) + (hn[, 2] * H[, 2]) + (hn[, 3] * H[, 3])
}



# gridHyper - calculate a grid for contouring.
#   Input:
#     rphi           = array (R, phi) ellipse axial ratios
#     kappa          = weighting parameter
#     nnodes         = number of grid nodes, n, in x and y
#     opts.normalize = normalize by n
#   Output:
#     z              = matrix of z values at the nxn grid nodes.
gridHyper <- function(rphi, rmax, kappa, nnodes, normalize = TRUE, proj = "eqd") {
  n <- nrow(rphi)
  stopifnot(n >= 2L)
  # Pre-allocate matrix (fixes subscript OOB)
  z <- matrix(0, nrow = nnodes, ncol = nnodes)

  s <- .R2zeta(rmax, proj)
  dx <- (2 * s) / (nnodes - 1)
  dy <- dx
  f <- if (normalize) kappa / (n^(1 / 3)) else kappa

  H <- .Rphi2H(rphi[, 1], rphi[, 2])

  # Precompute grid coordinates (x_i, y_j)
  xs <- seq(-s, s, length.out = nnodes)
  ys <- seq(-s, s, length.out = nnodes)

  # vectorized across data points using matrix operations for inner product
  for (i in seq_len(nnodes)) {
    xi <- xs[i]
    # build all y values for column i
    for (j in seq_len(nnodes)) {
      yj <- ys[j]
      # back-project node to (r,phi)
      rn_phi <- .xy2Rphi(xi, yj, rmax, proj)
      rn <- rn_phi[, "R"]
      pn <- rn_phi[, "phi"]
      hn <- .Rphi2H(rn, pn)
      hn_neg <- -hn
      # compute dot products with all H rows in a vectorized way:
      dvec <- .dotH(hn_neg, H) # length n
      zsum <- sum(exp(f * (1 - dvec)))
      z[i, j] <- zsum
    }
  }
  z
}

.drawFrame <- function(frame, proj, n = 512) {
  if (!frame) {
    return(NULL)
  }
  if (proj == "rfp") {
    data.frame(
      x = c(-1, 1, 1, 1, 1, -1, -1, -1),
      y = c(-1, -1, -1, 1, 1, 1, 1, -1)
    )
  } else {
    t <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cos(t), y = sin(t))
  }
}



#' R/phi plot
#'
#' @inheritParams mean_strain_ellipse
#' @inheritParams hypercontour
#' @param ... optional arguments passed to [hypercontour()]
#' @param contour logical. Whether a grid of densities should be drawn in the background.
#' @param image logical. Whether a raster image or filled contour lines should be
#' drawn as density grid. Ignored when `contour=FALSE`
#' @param n integer. Grid resolution or number of fille dcontours
#' @param contour.lines logical. Whether contour lines should be added.
#' @param contour.lines.params list of plotting arguments passed to [graphics::contour()]
#' @param col.palette function to produce color palette used for contouring
#' @param col.params list of plotting arguments passed `col.palette`
#' @param mean.ellipse logical. Whether the mean ellipse should be plotted
#' @param mean.ellipse.params list of plotting arguments passed to [ellipse()]
#' @param point.params list of plotting arguments passed to [graphics::points()]
#' @param at.x,at.y the points at which tick-marks and labels for the x and y
#' axes are to be drawn.
#'
#' @returns plot
#'
#' @seealso [hypercontour()], [Rphi_polar_plot()]
#'
#' @export
#'
#' @examples
#' data(ramsay)
#' Rphi_plot(ramsay[, 1], ramsay[, 2])
Rphi_plot <- function(r, phi,
                      contour = TRUE,
                      image = FALSE,
                      n = 10L,
                      contour.lines = TRUE,
                      contour.lines.params = list(lwd = .5, col = par("col")),
                      contour.col = viridis,
                      contour.col.params = list(),
                      mean.ellipse = TRUE,
                      mean.ellipse.params = list(border = "red", lwd = 2),
                      point.params = list(col = "grey", pch = 16, cex = .5),
                      rmax = NULL,
                      at.x = seq(-90, 90, 30),
                      at.y = NULL,
                      ...) {
  if (is.null(rmax)) {
    rmax <- ceiling(max(r)) + 1
  }

  plot(
    x = c(-90, 90),
    y = c(1, rmax),
    type = "n",
    axes = FALSE,
    xlim = c(-90, 90),
    ylim = c(1, rmax),
    xlab = bquote("Long-axis orientation," ~ varphi * " (" * degree * ")"),
    ylab = bquote("Strain ratio," ~ R[f]),
    main = "Rf/phi plot"
  )

  if (contour | contour.lines) {
    out <- hypercontour(r = r, phi = phi, proj = "rfp", rmax = rmax, ...)
  }


  if (contour) {
    if (image) {
      contour.col.params <- append(list(n = n), contour.col.params)
      colpal <- do.call(contour.col, contour.col.params)

      graphics::image(out$x, out$y, out$z,
        col = colpal,
        add = TRUE
      )
    } else {
      levels <- pretty(range(out$z, na.rm = TRUE), n)
      contour.col.params <- append(list(n = length(levels) - 1), contour.col.params)
      colpal <- do.call(contour.col, contour.col.params)

      graphics::.filled.contour(out$x, out$y, out$z, levels = levels, col = colpal)
    }
  }

  if (is.null(at.y)) at.y <- pretty(c(1, rmax))

  graphics::axis(1, at = at.x)
  graphics::axis(2, at = at.y)

  if (contour.lines) {
    contour.lines.params <- append(list(
      out$x, out$y, out$z,
      add = TRUE
    ), contour.lines.params)

    do.call(graphics::contour, contour.lines.params)
  }

  point.params <- append(list(x = out$points[, 1], y = out$points[, 2]), point.params)
  do.call(graphics::points, point.params)

  # add mean strain ellipse and its #95% confidence interval
  if (mean.ellipse) {
    rphi_mean <- mean_strain_ellipse(r, phi)

    mean.ellipse.params <- append(
      list(
        x = rphi_mean$phi, y = rphi_mean$R,
        radius.x = diff(rphi_mean$phi_CI) / 2,
        radius.y = diff(rphi_mean$R_CI) / 2
      ),
      mean.ellipse.params
    )

    do.call(ellipse, mean.ellipse.params)

    graphics::title(sub = bquote("R" == .(round(rphi_mean$R, 2)) ~ "|" ~ phi == .(round(rphi_mean$phi, 2)) * degree))
  }
}


#' Polar R/phi plot
#'
#' Plots R/phi fabric plot in polar coordinates
#'
#' @inheritParams mean_strain_ellipse
#' @param proj character. Projection,  `'eqd'` for equidistant (Elliot plot),
#' `'eqa'` for equal-area, `'stg'` for stereographic, `'ort'` for orthographic,
#' `'gno'` for gnomonic, `'lin'` for exponential (linear R), `'rdl'` for radial.
#' @param proj character.
#' @inheritParams Rphi_plot
#' @returns plot
#'
#' @seealso [hypercontour()], [Rphi_plot()]
#'
#' @export
#'
#' @examples
#' Rphi_polar_plot(ramsay[, 1], ramsay[, 2], proj = "eqa")
Rphi_polar_plot <- function(r, phi,
                            proj = c("eqd", "eqa", "stg", "ort", "gno", "lin", "rdl"),
                            contour = TRUE,
                            image = FALSE,
                            n = 10L,
                            contour.lines = TRUE,
                            contour.lines.params = list(lwd = .5, col = par("col")),
                            contour.col = viridis,
                            contour.col.params = list(),
                            mean.ellipse = TRUE,
                            mean.ellipse.params = list(col = "red", lwd = 2),
                            point.params = list(col = "grey", pch = 16, cex = .5),
                            ...) {
  proj <- match.arg(proj)

  out <- hypercontour(r, phi, proj = proj, ...)
  rmax <- max(out$rmax)

  plot(out$x, y = out$y, "n", xlab = NULL, ylab = NULL, asp = 1, axes = FALSE, ann = FALSE)

  if (contour) {
    if (image) {
      contour.col.params <- append(list(n = n), contour.col.params)
      colpal <- do.call(contour.col, contour.col.params)
      graphics::image(out$x, out$y, out$z,
        col = cols,
        add = TRUE
      )
    } else {
      levels <- pretty(range(out$z, na.rm = TRUE), n)
      contour.col.params <- append(list(n = length(levels) - 1), contour.col.params)
      colpal <- do.call(contour.col, contour.col.params)
      graphics::.filled.contour(out$x, out$y, out$z, levels = levels, col = colpal)
    }
  }

  if (contour.lines) {
    contour.lines.params <- append(list(
      out$x, out$y, out$z,
      add = TRUE
    ), contour.lines.params)

    do.call(graphics::contour, contour.lines.params)
  }

  point.params <- append(list(x = out$points[, 1], y = out$points[, 2]), point.params)
  do.call(graphics::points, point.params)

  if (mean.ellipse) {
    rphi_mean <- mean_strain_ellipse(r, phi)

    rphi_mean_xy <- .Rphi2xy(rphi_mean$R, deg2rad(rphi_mean$phi), rmax = rmax, proj = proj)

    mean.ellipse.params <- append(
      list(
        x = rphi_mean_xy[, 1],
        y = rphi_mean_xy[, 2]
      ),
      mean.ellipse.params
    )

    do.call(graphics::points, mean.ellipse.params)
    graphics::title(sub = bquote("R" == .(round(rphi_mean$R, 2)) ~ "|" ~ phi == .(round(rphi_mean$phi, 2)) * degree))
  }
  graphics::lines(out$frame, lwd = 2)



  proj.pretty <- c(
    "eqd" = "Equidistant",
    "eqa" = "Equal-area",
    "stg" = "Stereographic",
    "ort" = "Ortographic",
    "gno" = "Gnomonic",
    "lin" = "Linear",
    "rdl" = "Radial"
  )[proj]
  graphics::title(main = "Polar R/phi plot")
  graphics::mtext(paste(proj.pretty, "projection"))
}




mean_vorticity <- function(R) {
  (R^2 - 1) / (R^2 + 1)
}

crit_angle <- function(w, b = w) {
  (1 / 2) * asind(w / b) * (sqrt(1 - w^2) - sqrt(b^2 - w^2))
}

.b2r <- function(b) {
  sqrt((-b - 1) / (b - 1))
}

shape_factor <- function(r) {
  mean_vorticity(r)
}

.near <- function(x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}


RGN_hyperbola <- function(steps = 0.05, w = seq(0.1, 1, steps / 100)) {
  hyperbola_crit <- data.frame(
    b = w,
    theta = crit_angle(w),
    r = .b2r(w)
  )

  hyperbola_in <- expand.grid(w = w, b = w)
  hyperbola_in$theta <- crit_angle(hyperbola_in$w, hyperbola_in$b)
  hyperbola_in$r <- .b2r(hyperbola_in$b)

  hyperbola_in <- subset(hyperbola_in, !is.nan(hyperbola_in$theta) &
    (.near(hyperbola_in$w %% steps, 0) | .near(hyperbola_in$w %% steps, steps)))

  res <- split(hyperbola_in, hyperbola_in$w) |>
    lapply(function(h) {
      hs <- subset(h, theta == max(theta) &
        b == min(b))
      hs$theta <- 90
      hs
    })

  res2 <- rbind(hyperbola_in, do.call(rbind, res))
  hyperbola <- res2[order(res2$w, res2$b, -res2$theta), ]

  list(hyperbola = hyperbola, crit = hyperbola_crit)
}


vorticity_boot <- function(B, R = 100, probs = 0.975) {
  vapply(1:R, function(r) {
    quantile(sample(B, replace = TRUE), probs = probs, na.rm = TRUE) # take the upper 97.5% quantile to remove outliers
  }, FUN.VALUE = numeric(1))
}


#' Rigid Grain Net (RGN)
#'
#' The rigid grain net after (Jessup et al. 2007) plots the distribution the
#' strain ratio (`R`) of orientation (`phi`) of
#' porphyroclast over the theoretical distribution of tailless clasts. The plot estimates
#' the critical shape factor `Rc` marking the transition between the stable-sink
#' position and infinitely rotating porphyroclasts.
#' This critical shape factor can be interpreted as the the **mean kinmatic vorticity number**.
#' Here the `Rc` is estimated using the bootstrap method described in Stephan et al. (2025).
#'
#' @param r numeric. The porphyroclast aspect ratio (long axis/short axis)
#' @param theta numeric. Angle between long axis and foliation (in degrees)
#' @param angle_error numeric. Uncertainty of angle measurement. `3` by default.
#' @param probs integer. Probability with values in \eqn{[0, 1]} to estimate
#' critical shape factor, i.e. the largest shape factor of measurements outside
#' of critical hyperbole.
#' @param boot integer. Number of bootstrap resamples
#' @param grid numeric. Spacing of hyperboles.
#' @param ... plotting arguments passed to [graphics::points()]
#'
#' @returns a plot or a list of the calculated `B` (shape factor) and `theta` values,
#' and the bootstrapped confidence interval of the critical B value (`Rc_CI`).
#'
#' @references Jessup, Micah J., Richard D. Law, and Chiara Frassi.
#' "The rigid grain net (RGN): an alternative method for estimating mean
#' kinematic vorticity number (Wm)." Journal of Structural Geology 29.3 (2007):
#' 411-421. doi: 10.1016/j.jsg.2006.11.003
#'
#' Stephan, Tobias, et al. "Going with the flow—Changes of vorticity control
#' gold enrichment in Archean shear zones (Shebandowan Greenstone Belt,
#' Superior Province, Canada)." Journal of Structural Geology (2025): 105542.
#' doi: 10.1016/j.jsg.2025.105542
#'
#' @export
#'
#' @examples
#' data(ramsay)
#' set.seed(20250411)
#'
#' # assuming the mean orientation is the foliation:
#' theta <- tectonicr::circular_mean(ramsay[, 2]) - ramsay[, 2]
#'
#' RGN_plot(ramsay[, "R"], theta, col = "darkred")
RGN_plot <- function(r, theta, angle_error = 3, boot = 100L, probs = 0.972, grid = 0.05, ...) {
  R_val <- r
  theta <- theta %% 180
  theta <- ifelse(theta > 90, theta - 180, theta)

  B_val <- (R_val^2 - 1) / (R_val^2 + 1)
  # e_val = log(R_val) / 2


  crit_angleB <- crit_angle(B_val)
  infinite_rot_pos <- theta - angle_error > crit_angleB
  infinite_rot_neg <- theta + angle_error < crit_angleB
  crit <- abs(theta) - angle_error > crit_angleB


  bmax_r <- vorticity_boot(B_val[crit], R = boot, probs = probs)
  bmax_log <- log(bmax_r)

  bmax_geomean <- exp(mean(bmax_log, na.rm = TRUE))
  bmax_geosd <- exp(sd(bmax_log, na.rm = TRUE))


  R_test <- 10000
  t_score <- qt(p = 0.05 / 2, df = R_test - 1, lower.tail = FALSE)
  geo.sde <- bmax_geosd / sqrt(R_test)
  geo.margin_error <- t_score * geo.sde
  geo.lowerCI <- bmax_geomean - geo.margin_error
  geo.upperCI <- bmax_geomean + geo.margin_error

  hyp <- RGN_hyperbola(steps = grid)
  hyperbola <- split(hyp$hyperbola, hyp$hyperbola$w)

  plot(c(0, max(hyp$crit$b)), c(-90, 90),
    type = "n",
    axes = FALSE, frame.plot = FALSE,
    # xgap.axis = 0, ygap.axis = 0,
    ylim = c(-90, 90),
    xlab = "Shape factor, B*",
    ylab = expression("Angle between clast long axis and foliation," ~ theta ~ "(" * degree * ")")
  )

  # CI of critical B
  graphics::rect(geo.lowerCI, -100, geo.upperCI, 100, col = "grey80", border = NA)

  graphics::axis(2, at = seq(-90, 90, 30), gap.axis = 0)
  graphics::axis(1, at = seq(0, max(hyp$crit$b), .1), gap.axis = 0)


  # add hyperbola net
  lapply(hyperbola, function(h) {
    graphics::lines(h$b, h$theta, col = "grey85")
    graphics::lines(h$b, -h$theta, col = "grey85")
  })
  graphics::lines(hyp$crit$b, hyp$crit$theta, col = "grey30")
  graphics::lines(hyp$crit$b, -hyp$crit$theta, col = "grey30")

  graphics::abline(v = cosd(45), col = "grey10", lty = 3, lwd = .1) # pure-shear simple shear separation
  graphics::abline(v = c(geo.lowerCI, geo.upperCI), col = "black", lty = 2, lwd = .5) # CI interval of bootstrapped B

  graphics::points(B_val, theta, ...)


  graphics::mtext(paste0("Rc = ", round(bmax_geomean, 2), " \u00B1 ", round(geo.margin_error, 2)))
  graphics::title(sub = paste0("(n: ", length(r), ")"))


  invisible(
    list(
      values = cbind(B = B_val, theta = theta),
      Rc_CI = c(geo.lowerCI, geo.upperCI)
    )
  )
}
