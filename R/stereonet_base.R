#' @keywords internal
roty3 <- function(deg) {
  theta <- deg2rad(deg)
  c <- cos(theta)
  s <- sin(theta)

  matrix(c(
    c, 0, -s,
    0, 1,  0,
    s, 0,  c
  ), nrow = 3, byrow = TRUE)
}

#' @keywords internal
rotz3 <- function(deg) {
  theta <- deg2rad(deg)
  c <- cos(theta)
  s <- sin(theta)

  matrix(c(
    c, s, 0,
    -s, c, 0,
    0, 0, 1
  ), nrow = 3, byrow = TRUE)
}


#' @keywords internal
.fix_inc <- function(az, inc) {
  tinc <- deg2rad(inc %% 360)
  co <- cos(tinc)
  si <- sin(tinc)
  ang <- atan2d(si, co)
  quad <- rep(1, length(inc))
  quad[co >= 0 & si >= 0] <- 1
  quad[co < 0 & si >= 0] <- 2
  quad[co < 0 & si < 0] <- 3
  quad[co >= 0 & si < 0] <- 4
  tinc[quad == 1] <- ang[quad == 1]
  tinc[quad == 2] <- 180 - ang[quad == 2]
  tinc[quad == 3] <- 180 + ang[quad == 3]
  tinc[quad == 4] <- -ang[quad == 4]
  az[quad == 1] <- az[quad == 1]
  az[quad == 2] <- 180 + az[quad == 2]
  az[quad == 3] <- az[quad == 3]
  az[quad == 4] <- 180 + az[quad == 4]
  az <- az %% 360
  inc <- tinc
  return(cbind(az, inc))
}




#' Stereographic projection
#'
#' Transformation of spherical coordinates into the stereographic projection
#'
#' @param az,inc numeric vectors. Azimuth and Inclination in degrees.
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param r numeric. Radius of circle. Default is `1` for unit circle.
#' @param earea logical `TRUE` for Lambert equal-area projection (also "Schmidt net"; the default), or
#' `FALSE` for meridional stereographic projection (also "Wulff net" or "Stereonet").
#'
#' @returns two-column vector with the transformed coordinates
#'
#' @export
#'
#' @examples
#' stereo_coords(90, 10)
#' stereo_coords(90, 10, earea = TRUE, upper.hem = TRUE)
stereo_coords <- function(az, inc, upper.hem = FALSE, earea = TRUE, r = 1) {
  if (upper.hem) {
    az <- az + 180
  }


  B <- .fix_inc(az = az, inc = 90 - inc)
  azi <- deg2rad(B[, 1])
  inc <- deg2rad(B[, 2])

  if (earea) {
    tq <- sqrt(2) * sin(inc / 2)
  } else {
    tq <- tan(inc / 2)
  }
  pltx <- r * tq * sin(azi)
  plty <- r * tq * cos(azi)
  cbind(x = pltx, y = plty)
}

#' @keywords internal
.schmidt_crds <- function(az_rad, inc_rad, r = 1) {
  tq <- r * sqrt(2) * sin(pi / 4 - inc_rad / 2)
  x <- tq * sin(az_rad)
  y <- tq * cos(az_rad)
  cbind(x = x, y = y)
}

#' @keywords internal
.wulff_crds <- function(az_rad, inc_rad, r = 1) {
  tq <- r * tan(pi / 4 - inc_rad / 2)
  x <- tq * sin(az_rad)
  y <- tq * cos(az_rad)
  cbind(x = x, y = y)
}


#' Stereographic projection of lines and planes
#'
#' Visualization of lines, planes in a stereographic projection.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`. `"Pair"`, or `"Fault"`
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param earea logical `TRUE` for Lambert equal-area projection (also "Schmidt net"; the default), or
#' `FALSE` for meridional stereographic projection (also "Wulff net" or "Stereonet").
#' @param col color
#' @param pch plotting character
#' @param lab character. text labels
#' @param text.pos position for labels
#' @param cex character expansion of labels
#' @param ... optional graphical parameters
#' @note `"Plane"` and `"Fault"` objects will be displayed as pole to the plane (only).
#' @importFrom graphics points text
#' @export
#' @examples
#' stereoplot()
#' stereo_point(Line(c(90, 80), c(10, 75)), lab = c("L1", "L2"))
#' stereo_point(Plane(120, 30), lab = "P", col = "red")
stereo_point <- function(x, col = 1, pch = 20, lab = NULL, text.pos = 4, cex = 1, upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.spherical(x))
  if (is.Vec3(x)) x <- Line(x)

  if (is.Plane(x) | is.Fault(x)) {
    x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  }

  crds <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem, earea
  )

  graphics::points(crds[, "x"], crds[, "y"], pch = pch, col = col, cex = cex, ...)
  if (!is.null(lab)) {
    graphics::text(crds[, "x"], crds[, "y"], labels = lab, pos = text.pos, col = col)
  }
}


#' Stereographic projection of pairs
#'
#' Visualization of pairs (planes and lines) in a stereographic projection.
#'
#' @param x Object of class `"Fault"`
#' @param greatcircles logical. Whether greatcircles are displayed (`TRUE`, the default) or poles to planes (`FALSE`)
#' @param pch,col,lwd,lty plotting parameters for planes and lines
#' @param lab character. text labels
#' @param cex character expansion of labels
# #' @param text.pos position for labels
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param earea logical `TRUE` for Lambert equal-area projection (also "Schmidt net"; the default), or
#' `FALSE` for meridional stereographic projection (also "Wulff net" or "Stereonet").
#'
#' @note `"Plane"` objects will be displayed as pole to the plane.
#' @importFrom graphics points text
#' @name stereo-pair
#'
#' @examples
#' faults <- Fault(
#'   c(0, 90, 180, 270),
#'   c(80, 45, 80, 45),
#'   c(0, 170, 180, 315),
#'   c(80, 10, 80, 36),
#'   c(1, -1, 0, 1)
#' )
#' # stereoplot()
#' # stereo_fault(faults, col = 1:4)
#' # stereo_fault(faults, col =1:4, hoeppener = TRUE)
#' # legend("bottomright", c("normal", "thrust", "unknown", "normal"), fill = 1:4)
#'
#' stereoplot()
#' stereo_pair(faults, col = 1:4)
NULL

# #' @rdname stereo-fault
# #' @export
# stereo_fault <- function(x, hoeppener = FALSE, greatcirles = TRUE, pch = 16, col = 1, lwd = 1, lty = 1, lab = NULL, cex = 1, text.pos = 4, upper.hem = FALSE, earea = TRUE) {
#   stopifnot(is.Fault(x))
#   x0 <- x
#
#   if (length(col) == 1) col <- rep(col, nrow(x))
#   if (length(pch) == 1) pch <- rep(pch, nrow(x))
#   if (length(cex) == 1) cex <- rep(cex, nrow(x))
#   if (length(lwd) == 1) lwd <- rep(lwd, nrow(x))
#   if (length(lty) == 1) lty <- rep(lty, nrow(x))
#
#   x[, 1] <- 180 + x[, 1]
#   x[, 2] <- 90 - x[, 2]
#
#   crds.p <- stereo_coords(
#     x[, 1],
#     x[, 2],
#     upper.hem
#   )
#   if (hoeppener) {
#     crds.l <- stereo_coords(
#       x[, 3] + 180,
#       90 - x[, 4],
#       upper.hem, earea
#     )
#   } else {
#     crds.l <- stereo_coords(
#       x[, 3],
#       x[, 4],
#       upper.hem, earea
#     )
#   }
#
#   for (i in 1:nrow(crds.p)) {
#     if (greatcirles) {
#       stereo_greatcircle(Plane(x0[i, 1], x0[i, 2]), col = col[i], lwd = lwd[i], lty = lty[i])
#     } else {
#       graphics::points(crds.p[i, "x"], crds.p[i, "y"], pch = pch[i], col = col[i])
#       if (!is.null(lab)) {
#         graphics::text(crds.p[, "x"], crds.p[, "y"], labels = lab, pos = text.pos, col = col[i], cex = cex[i])
#       }
#     }
#     if (hoeppener) {
#       ang <- (x0[i, 3] * sign(x0[i, "sense"]))
#       graphics::points(crds.l[i, "x"], crds.l[i, "y"], pch = pch[i], col = col[i], cex = cex[i] / 1.25)
#       if (x0[i, "sense"] != 0) {
#         graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = "\u2191", col = col[i], srt = ang, cex = cex[i] * 1.5)
#       }
#     } else {
#       ang <- x0[i, 3] %% 360
#
#       symb <- if (x0[i, 5] < 0) {
#         "\u2193"
#       } else {
#         "\u2191"
#       }
#
#       if (x0[i, "sense"] != 0) {
#         graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = symb, col = col[i], srt = -ang, cex = cex[i] * 1.5)
#       } else {
#         graphics::points(crds.l[i, "x"], crds.l[i, "y"], pch = pch[i], col = col[i], cex = cex[i])
#       }
#     }
#   }
# }

#' @rdname stereo-pair
#' @export
stereo_pair <- function(x, pch = 16, col = 1, lwd = 1, lty = 1, lab = NULL, cex = 1, greatcircles = TRUE, upper.hem = FALSE, earea = TRUE) {
  if (greatcircles) {
    stereo_greatcircle(Fault_plane(x), lwd = lwd, lty = lty, col = col, upper.hem = upper.hem, earea = earea)
  } else {
    stereo_point(Fault_plane(x), pch = pch, cex = cex, col = col, upper.hem = upper.hem, earea = earea)
  }

  stereo_point(Fault_slip(x), pch = pch, cex = cex, col = col, upper.hem = upper.hem, earea = earea)
}


#' Stereographic projection of cones
#'
#' Visualization of smallcircles and greatcircles in a stereographic projection.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @param d numeric. conical angle in degrees.
#' @param col,lty,lwd color, line type, and line width parameters
#' @param N integer. number of points to calculate
#' @param BALL.radius numeric size of sphere
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param earea logical `TRUE` for Lambert equal-area projection (also "Schmidt net"; the default), or
#' `FALSE` for meridional stereographic projection (also "Wulff net" or "Stereonet").
#' @param ... optional graphical parameters
#' @importFrom graphics lines
#' @name stereo_cones
#' @examples
#' stereoplot()
#' stereo_point(Line(90, 5), lab = "L")
#' stereo_smallcircle(Line(90, 5), d = 10)
#' stereo_point(Plane(120, 30), lab = "P", col = "red")
#' stereo_greatcircle(Plane(120, 30), col = "red")
#'
#' stereoplot()
#' stereo_point(Line(c(129, 90), c(30, 5)), lab = c("L1", "L2"))
#' stereo_smallcircle(Line(c(129, 90), c(30, 5)), d = c(10, 5), col = 1:2, lty = 1:2, lwd = 1:2)
NULL

#' @rdname stereo_cones
#' @export
stereo_smallcircle <- function(x, d = 90, col = 1, N = 1000, upper.hem = FALSE, earea = TRUE, lty = 1, lwd = 1, BALL.radius = 1, ...) {
  if (length(d) == 1) {
    stereo_smallcircle0(x, d, col, N, upper.hem, earea, lty, lwd, BALL.radius, ...)
  } else {
    nx <- nrow(x)
    stopifnot(length(d) == nx)

    if (length(col) == 1) col <- rep(col, nx)
    if (length(lwd) == 1) lwd <- rep(lwd, nx)
    if (length(lty) == 1) lty <- rep(lty, nx)

    invisible(
      lapply(seq_len(nx), function(i) {
        stereo_smallcircle0(Line(x[i, ]), d[i], col[i], N, upper.hem, earea, lty[i], lwd[i], BALL.radius, ...)
      })
    )
  }
}

stereo_smallcircle0 <- function(x, d = 90, col = 1, N = 1000, upper.hem = FALSE, earea = TRUE, lty = 1, lwd = 1, BALL.radius = 1, ...) {
  stopifnot(is.spherical(x))
  if (length(col) == 1) col <- rep(col, nrow(x))
  if (length(lty) == 1) lty <- rep(lty, nrow(x))
  if (length(lwd) == 1) lwd <- rep(lwd, nrow(x))

  if (is.Vec3(x)) x <- Line(x)
  az <- x[, 1]
  inc <- 90 - x[, 2]

  phi <- seq(from = 0, to = 2 * pi, length = N)
  theta <- deg2rad(d)

  x <- BALL.radius * sin(theta) * cos(phi)
  y <- BALL.radius * sin(theta) * sin(phi)
  z <- rep(BALL.radius * cos(theta), N)
  D <- cbind(x, y, z)

  for (i in 1:length(az)) {
    ry <- roty3(inc[i])
    rz <- rotz3(az[i])
    Rmat <- ry %*% rz
    g <- D %*% Rmat
    r2 <- sqrt(g[, 1]^2 + g[, 2]^2 + g[, 3]^2)
    phi2 <- atan2d(g[, 2], g[, 1])
    theta2 <- acosd(g[, 3] / r2)
    Sc <- stereo_coords(phi2, 90 - theta2, upper.hem, earea)

    diss <- sqrt((Sc[1:(N - 1), "x"] - Sc[2:(N), "x"])^2 + (Sc[1:(N - 1), "y"] - Sc[2:(N), "y"])^2)
    ww <- which(diss > 0.9 * BALL.radius)
    if (length(ww) > 0) {
      Sc[ww, "x"] <- NA
      Sc[ww, "y"] <- NA
    }
    graphics::lines(Sc[, "x"], Sc[, "y"], col = col[i], lty = lty[i], lwd = lwd[i], ...)
  }
}

#' @rdname stereo_cones
#' @export
stereo_greatcircle <- function(x, ...) {
  stopifnot(is.spherical(x))
  x[, 2] <- 90 + x[, 2]
  stereo_smallcircle(x, d = 90, ...) # add circle
}


#' Great-circle segment between two vectors
#'
#' Plots the great-circle segment between two vectors
#'
#' @param x,y objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`
#' @inheritParams stereo_smallcircle
#' @param n integer. number of points along greatcircle (100 by default)
#' @param ... graphical parameters passed to [graphics::lines()]
#'
#' @returns NULL
#' @export
#'
#' @seealso [slerp()], [stereo_greatcircle]
#'
#' @examples
#' x <- Line(120, 7)
#' y <- Line(10, 13)
#' plot(rbind(x, y))
#' stereo_segment(x, y, col = "red")
#'
#' # For multiple segments use lapply():
#' set.seed(20250411)
#' mu <- Line(45, 10)
#' x <- rvmf(100, mu = mu)
#' plot(x)
#' lapply(seq_len(nrow(x)), FUN = function(i) {
#'   stereo_segment(x[i, ], mu, col = i)
#' })
#' points(mu, pch = 16, col = "white")
stereo_segment <- function(x, y, upper.hem = FALSE, earea = TRUE, n = 100L, BALL.radius = 1, ...) {
  stopifnot(nrow(x) == 1, nrow(y) == 1)
  vang <- angle(x, y)
  if (is.Vec3(x)) vang <- rad2deg(vang)

  if (vang <= 90) {
    n_part <- vang / 90
    n <- ceiling(n * n_part)

    .draw_lines(x, y, n, upper.hem, earea, BALL.radius, ...)
  } else {
    xy <- crossprod(x, y) |> Line()
    strike1 <- Line(xy[1, 1] + 90, 0)
    strike2 <- Line(xy[1, 1] - 90, 0)

    n_part <- (180 - vang) / 180

    na <- ceiling(n * n_part)
    nb <- ceiling(n * (1 - n_part))

    if (angle(x, strike1) <= angle(x, strike2)) {
      line1 <- strike1
      line2 <- strike2
      n1 <- nb
      n2 <- na
    } else {
      line1 <- strike2
      line2 <- strike1
      n1 <- na
      n2 <- nb
    }

    .draw_lines(x, line1, n1, upper.hem, earea, BALL.radius, ...)
    .draw_lines(line2, y, n2, upper.hem, earea, BALL.radius, ...)
  }
}

.draw_lines <- function(x, y, n = 100L, upper.hem, earea, BALL.radius = 1, ...) {
  t <- seq(0, 1, length.out = n)
  D <- slerp(x, y, t) |>
    Line() |>
    unclass()

  Sc <- stereo_coords(D[, 1], D[, 2], upper.hem, earea)

  diss <- sqrt((Sc[1:(n - 1), "x"] - Sc[2:(n), "x"])^2 + (Sc[1:(n - 1), "y"] - Sc[2:(n), "y"])^2)
  ww <- which(diss > 0.9 * BALL.radius)
  if (length(ww) > 0) {
    Sc[ww, "x"] <- NA
    Sc[ww, "y"] <- NA
  }
  graphics::lines(Sc[, "x"], Sc[, "y"], ...)
}


#' Stereoplot frame
#'
#' Adds a (primitive) circle with given radius to an existing plot
#'
#' @param n integer. Resolution of circle's line
#' @inheritParams stereoplot
#' @param ... optional arguments passed to [graphics::lines()]
#'
#' @export
#' @seealso [stereoplot()], [stereoplot_ticks()], [stereoplot_guides()]
#'
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_frame(col = "red", lwd = 3)
stereoplot_frame <- function(n = 512L, radius = 1, ...) {
  phi <- seq(0, 2 * pi, by = 2 * pi / n)
  x <- cos(phi) * radius
  y <- sin(phi) * radius
  graphics::lines(x, y, ...)
}


#' Stereographic projection
#'
#' Initialize the plot for equal-area stereographic projections (Wulff) or Lambert
#' Equal-Area projections (Schmidt).
#'
#' @param earea logical. Projection, either `TRUE` for Lambert equal-area
#' projection (the default), or `FALSE` for meridional stereographic projection.
#' @param guides logical. Whether guides should be added to the plot (`TRUE` by default)
#' @param d integer. Angle distance between guides. Default: 10
#' @param col color of guide lines
#' @param lwd linewidth of guide lines
#' @param lty linetype of guide lines
#' @param border.col color of primitive circle (frame), center-cross and ticks of the stereo plot
#' @param centercross logical. Whether a center cross should be added (`TRUE` by default)
#' @param ticks integer. Angle between ticks. if `NULL` (the default), no ticks are drawn.
#' @param title,sub character. Title and subtitle of plot
#' @param origin.text character. Text at origin of stereoplot.
#' @param labels this can either be a logical value specifying whether (numerical)
#' annotations are to be made next to the tickmarks, or a character or expression
#' vector of labels to be placed next to the tickpoints.
#' @param ladj adjustment for all labels away from origin of stereoplot circle.
#' This essentially an amount that is added to `radius` and the length of the ticks.
#' @param radius numeric. Radius of circle
#'
#' @source Adapted from the `RFOC` package
#' @importFrom graphics points
#' @export
#' @examples
#' stereoplot(ticks = 30, title = "title", sub = "subtitle", border.col = "purple", labels = TRUE)
stereoplot <- function(earea = TRUE, guides = TRUE, d = 10, col = "lightgray",
                       lwd = 1, lty = 1, border.col = "black", title = NULL,
                       sub = NULL, origin.text = "N", labels = FALSE, ladj = 0.04, centercross = TRUE, ticks = NULL, radius = 1) {
  plot(radius * c(-1, 1), radius * c(-1, 1),
    type = "n", xlab = NULL, ylab = NULL, asp = 1,
    axes = FALSE, ann = FALSE
  )

  graphics::title(main = title, sub = sub)
  graphics::mtext(origin.text, col = border.col)

  if (guides) stereoplot_guides(d = d, earea = earea, col = col, lwd = lwd, lty = lty, radius = radius)

  if (!is.null(ticks)) stereoplot_ticks(angle = ticks, col = border.col, radius = radius, labels = labels, ladj = ladj)

  if (centercross) points(0, 0, pch = 3, col = border.col)

  stereoplot_frame(col = border.col, radius = radius)
}

#' Stereoplot tickmarks
#'
#' Adds stereoplot rickmarks to an existing plot
#'
#' @inheritParams stereoplot
#' @param length numeric. Length of ticks as fraction of `radius`
#' @param angle numeric. Division angle in degrees
#' @param rotation numeric. Rotation (positive for counter-clockwise) of tickmarks and labels
#' @param ... optional arguments passed to [graphics::segments()] and [graphics::text()]
#'
#' @importFrom graphics segments
#' @export
#' @seealso [stereoplot()], [stereoplot_frame()], [stereoplot_guides()]
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_frame()
#' stereoplot_ticks(length = 0.05, angle = 45, col = "blue", lwd = 2, labels = TRUE)
stereoplot_ticks <- function(length = 0.02, angle = 10, labels = FALSE, ladj = 2 * length, radius = 1, rotation = 0, ...) {
  DR <- radius + length
  ang_deg <- (seq(0, 360, by = angle) + rotation) %% 360
  ang <- pi * ang_deg / 180
  graphics::segments(
    radius * cos(ang), radius * sin(ang), DR * cos(ang), DR * sin(ang),
    ...
  )

  do.label <- is.logical(labels)
  if (do.label & isTRUE(labels)) {
    labels <- seq(0, 360, by = angle)
    labels[length(labels)] <- NA
    labels <- as.character((-labels + 90) %% 360)
  } else {
    labels <- NULL
  }
  if (!is.null(labels)) {
    stopifnot(length(labels) == length(ang))
    DR.labs <- DR + ladj
    graphics::text(
      DR.labs * cos(ang), DR.labs * sin(ang),
      labels = labels, ...
    )
  }
}



stereo_guides_schmidt <- function(d = 10, n = 512, r = 1, rotation = 0, ...) {
  rot <- deg2rad(rotation)
  lam_seq <- deg2rad(seq(0, 180, length.out = n))
  lam0 <- pi / 2
  R <- sqrt(2) / 2

  # Precompute sin and cos of lam - lam0
  cos_lam <- cos(lam_seq - lam0)
  sin_lam <- sin(lam_seq - lam0)

  # Latitude lines (constant phi)
  phi_vals <- deg2rad(seq(-90 + d, 90 - d, 10))
  cos_phi <- cos(phi_vals)
  sin_phi <- sin(phi_vals)

  for (i in seq_along(phi_vals)) {
    phi <- phi_vals[i]
    kp <- sqrt(2 / (1 + cos_phi[i] * cos_lam))
    x <- R * kp * cos_phi[i] * sin_lam
    y <- R * kp * sin_phi[i]
    graphics::lines(x, y, ...)
  }

  # Longitude lines (constant lambda)
  phi_seq <- deg2rad(seq(-90, 90, 5))
  cos_phi_seq <- cos(phi_seq)
  sin_phi_seq <- sin(phi_seq)

  lam_vals <- deg2rad(seq(d, 180 - d, d))

  lam_vals_rot <- (lam_vals - lam0)
  cos_lam_vals <- cos(lam_vals_rot)
  sin_lam_vals <- sin(lam_vals_rot)

  for (j in seq_along(lam_vals)) {
    lam <- lam_vals[j]
    kp <- sqrt(2 / (1 + cos_phi_seq * cos_lam_vals[j]))
    x <- R * kp * cos_phi_seq * sin_lam_vals[j] * r
    y <- R * kp * sin_phi_seq * r
    graphics::lines(x, y, ...)
  }
}

stereo_guides_wulff <- function(d = 9, n = 512, r = 1, rotation = 0, ...) {
  rot <- deg2rad(rotation)
  if (n %% 2 != 0) n <- n + 1
  beta0 <- deg2rad(seq(0, 360, length.out = n))
  beta <- c(beta0[(n / 2):n], beta0[1:(n / 2 - 1)])

  cos_beta <- cos(beta)
  sin_beta <- sin(beta)

  # great circles
  phi <- deg2rad(seq(0, 180, by = d))
  rcos_phi <- r / cos(phi)
  rtan_phi <- r * tan(phi)


  for (i in seq_along(phi)) {
    xg <- -rtan_phi[i] + (rcos_phi[i]) * cos_beta
    yg <- rcos_phi[i] * sin_beta
    condg <- sqrt(xg^2 + yg^2) <= r
    graphics::lines(xg[condg], yg[condg], ...)
  }


  # small circles
  gamma <- phi
  rtan_gamma <- r * tan(gamma)
  cos_gamma <- cos(gamma)
  rcos_gamma <- r / cos_gamma
  for (j in seq_along(gamma)) {
    xs <- rtan_gamma[j] * cos_beta
    ys1 <- rcos_gamma[j] + rtan_gamma[j] * sin_beta
    ys2 <- -rcos_gamma[j] + rtan_gamma[j] * sin_beta

    conds1 <- sqrt(xs^2 + ys1^2) <= r
    conds2 <- sqrt(xs^2 + ys2^2) <= r

    graphics::lines(xs[conds1], ys1[conds1], ...)
    graphics::lines(xs[conds2], ys2[conds2], ...)
  }

  graphics::lines(c(0, 0), c(1, -1), ...)
  graphics::lines(c(-1, 1), c(0, 0), ...)
}

#' Stereoplot gridlines
#'
#' Adds equal-area or equal-angle projection gridlines to an existing stereoplot.
#'
#' @param d angle between grid lines
#' @inheritParams stereoplot
#' @param ... optional arguments passed to [graphics::lines()]
#'
#' @importFrom graphics lines
#'
#' @seealso [stereoplot()], [stereoplot_frame()], [stereoplot_ticks()]
#' @export
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_guides(d = 5, earea = FALSE, col = "green", rotation = 20)
#'
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_guides(d = 15, earea = TRUE, col = "orange", rotation = 90)
stereoplot_guides <- function(d = 10, earea = TRUE, radius = 1, ...) {
  if (earea) {
    stereo_guides_schmidt(d = d, r = radius, ...)
  } else {
    stereo_guides_wulff(d = d, r = radius, ...)
  }
}


#' Plot spherical objects
#'
#' @param x objects of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @inheritParams stereoplot
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param grid.params list.
#' @param ... parameters passed to [stereo_point()], [stereo_smallcircle()], [stereo_greatcircle()], or [fault_plot()]
#'
#' @details
#' If `x` is a Ray, than solid symbols show rays pointing in the lower hemisphere,
#' and open symbols point into the upper hemisphere
#'
#' @name plot-spherical
#'
#' @examples
#' plot(rvmf(10, mu = Vec3(1, 0, 0))) # Vec
#' plot(Line(c(90, 80), c(10, 75)), lab = c("L1", "L2"))
#' plot(Ray(c(90, 80), c(10, 75), sense = c(1, -1)), lab = c("L1", "L2"))
#' plot(Plane(120, 30), col = "red")
#' plot(Pair(120, 50, 36, 8))
#' plot(Fault(120, 50, 36, 8, -1))
NULL

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Line <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  stereo_point(x, upper.hem = upper.hem, earea = earea, ...)
}

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Vec3 <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  stereo_point(x, upper.hem = upper.hem, earea = earea, ...)
}

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Ray <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  if (is.Ray(x)) {
    # different symbols for upper and lower hemisphere rays
    pch_upp <- if (upper.hem) c(1, 16) else c(16, 1)
    stereo_point(x,
      upper.hem = upper.hem, earea = earea,
      pch = ifelse(x[, 2] > 1, pch_upp[1], pch_upp[2]),
      ...
    )
  }
}

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Plane <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  stereo_greatcircle(x, upper.hem = upper.hem, earea = earea, ...)
}

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Pair <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  stereo_greatcircle(Plane(x), upper.hem = upper.hem, earea = earea, ...)
  stereo_point(Line(x), upper.hem = upper.hem, earea = earea, ...)
}

#' @rdname plot-spherical
#' @exportS3Method graphics::plot
plot.Fault <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))
  fault_plot(x, upper.hem = upper.hem, earea = earea, ...)
}


# plot.spherical <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
#   do.call(stereoplot, append(grid.params, earea))
#
#   if (is.Line(x) | is.Vec3(x)) stereo_point(x, upper.hem = upper.hem, earea = earea, ...)
#   if(is.Ray(x)){
#     # different symbols for upper and lower hemisphere rays
#     pch_upp <- if(upper.hem) c(1, 16) else c(16, 1)
#     stereo_point(x, upper.hem = upper.hem, earea = earea,
#                  pch = ifelse(x[, 2] > 1, pch_upp[1], pch_upp[2]),
#                  ...)
#   }
#   if (is.Plane(x) & !is.Pair(x)) stereo_greatcircle(x, upper.hem = upper.hem, earea = earea, ...)
#   if (is.Fault(x)) fault_plot(x, upper.hem = upper.hem, earea = earea, ...)
#   if (is.Pair(x) & !is.Fault(x)) {
#     stereo_greatcircle(Fault_plane(x), upper.hem = upper.hem, earea = earea, ...)
#     stereo_point(Fault_slip(x), upper.hem = upper.hem, earea = earea, ...)
#   }
# }
#

#' Add Points to a Plot
#'
#' @param ... arguments passed to [graphics::points()]
#' @inheritParams plot.Vec3
#' @inheritParams graphics::text
#' @importFrom graphics points
#'
#' @exportS3Method graphics::points
#'
#' @examples
#' stereoplot()
#' points(rvmf(n = 100))
#' points(Plane(120, 30), col = "red", pch = 19)
points.spherical <- function(x, upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.spherical(x))
  if (is.Vec3(x)) x <- Line(x)

  if (is.Plane(x) | is.Fault(x)) {
    x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  }

  crds <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem, earea
  )

  graphics::points(crds[, "x"], crds[, "y"], ...)
}


#' Add Lines to a Plot
#'
#' @inheritParams plot.Vec3
#' @inheritParams graphics::text
#' @param ang numeric. Conical angle in degrees.
#' @param ... arguments passed to [graphics::lines()]
#' @importFrom graphics lines
#' @exportS3Method graphics::lines
#'
#' @examples
#' set.seed(20250411)
#' stereoplot()
#' lines(rvmf(n = 5), ang = runif(5, 0, 90), col = 1:5)
lines.spherical <- function(x, ang = 90, ...) {
  if (is.Plane(x)) stereo_greatcircle(x, ...) else stereo_smallcircle(x, d = ang, ...)
}


#' Add Points to a Plot
#'
#' @param ... arguments passed to [graphics::text()]
#' @inheritParams plot.Vec3
#' @inheritParams graphics::text
#' @importFrom graphics text
#'
#' @exportS3Method graphics::text
#'
#' @examples
#' stereoplot()
#' points(Line(c(90, 80), c(10, 75)), col = 1:2)
#' text(Line(c(90, 80), c(10, 75)), labels = c("L1", "L2"), col = 1:2, pos = 3)
text.spherical <- function(x, labels = seq_along(x[, 1]), upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.spherical(x))
  if (is.Vec3(x)) x <- Line(x)

  if (is.Plane(x) | is.Fault(x)) {
    x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  }

  crds <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem, earea
  )

  graphics::text(crds[, "x"], crds[, "y"], labels = labels, ...)
}


hypot <- function(x, y) {
  sqrt(x^2 + y^2)
}

#' Add Arrows to a Stereoplot
#'
#' A quiver plot displays displacement vectors into pointing into the direction of movement.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`.
#' @param sense numeric. Sense of the line on a fault plane. Either
#' `1`or `-1` for normal or thrust offset, respectively. The "sense" is the sign
#' of the fault's rake (see [Fault_from_rake()] for details).
#' @param ... arguments passed to [graphics::arrows()]
#' @inheritParams plot.Vec3
#' @param length numeric. Length of the edges of the arrow head (in inches).
#' @param angle numeric. Angle from the shaft of the arrow to the edge of the arrow head.
#' @param scale numeric. Scales the length of the vector. `0.1` by default
#'
#' @seealso [hoeppener()], [angelier()]
#'
#' @importFrom graphics arrows
#' @export
#'
#' @examples
#' set.seed(20250411)
#' stereoplot()
#' p <- rvmf(n = 100)
#' points(p, pch = 16, cex = .5)
#' stereo_arrows(p, sense = 1, col = "red")
stereo_arrows <- function(x, sense, scale = .1, angle = 10, length = 0.1, upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Ray(x) | is.Plane(x))

  if (nrow(x) > 1 & length(sense) == 1) sense <- rep(sense, nrow(x))

  if (is.Vec3(x)) x <- Line(x)

  if (is.Plane(x)) {
    x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  }

  crds <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem = upper.hem, earea
  )

  dx <- crds[, "x"]
  dy <- crds[, "y"]

  mag <- hypot(dx, dy)
  u <- sense * dx / mag
  v <- sense * dy / mag

  graphics::arrows(dx, dy, dx + scale * u, dy + scale * v, angle = angle, length = length, ...)
}


#' Add fault data to existing plot
#'
#' @param x `"Fault"` object
#' @param type character. One of `"angelier"` (for "Angelier plot") or `"hoeppener"` (for "Hoeppener plot"). See details.
#' @param lty,lwd,cex,pch,col,bg plotting parameters
#' @param points logical. Whether the lineation points (Angelier plot) or poles (Hoeppner plot) should be added to the plot
#' @param ... arguments passed to [stereo_arrows()]
#'
#' @returns NULL
#' @name fault-plot
#'
#' @seealso [stereo_arrows()]
#'
#' @details
#' **Angelier plot** shows all planes as *great circles* and lineations as points. Fault striae are plotted as vectors on top of the lineation pointing in the movement direction of the hangingwall. Easy to read in case of homogeneous or small datasets.
#'
#' **Hoeppener plot** shows all planes as *poles* while lineations are not shown. Fault striae are plotted as vectors on top of poles pointing in the movement direction of the hangingwall. Useful in case of large or heterogeneous datasets.
#'
#' @references
#' Angelier, J. Tectonic analysis of fault slip data sets, J. Geophys. Res. 89 (B7), 5835-5848 (1984)
#'
#' Hoeppener, R. Tektonik im Schiefergebirge. Geol Rundsch 44, 26-58 (1955). https://doi.org/10.1007/BF01802903
#'
#' @examples
#' f <- Fault(
#'   c("a" = 120, "b" = 125, "c" = 100),
#'   c(60, 62, 50),
#'   c(110, 25, 30),
#'   c(58, 9, 23),
#'   c(1, -1, 1)
#' )
#'
#' stereoplot(title = "Angelier plot")
#' angelier(f, col = 1:nrow(f), pch = 16, scale = 0.1)
#'
#' stereoplot(title = "Hoeppener plot")
#' hoeppener(f, col = 1:nrow(f), cex = 1, scale = 0.1, points = FALSE)
#'
#' # or
#' stereoplot()
#' fault_plot(f, type = "hoeppener", col = 1:nrow(f), cex = 1, scale = 0.1, points = FALSE)
NULL

#' @rdname fault-plot
#' @export
fault_plot <- function(x, type = c("angelier", "hoeppener"), ...) {
  type <- match.arg(type)
  if (type == "angelier") angelier(x, ...) else hoeppener(x, ...)
}

#' @rdname fault-plot
#' @export
hoeppener <- function(x, pch = 1, col = "black", cex = 1, bg = NULL, points = TRUE, ...) {
  stopifnot(is.Fault(x))

  stereo_arrows(Plane(x), sense = x[, "sense"], col = col, ...)
  if (isTRUE(points)) points(Plane(x), pch = pch, col = col, cex = cex, bg = bg)
}

#' @rdname fault-plot
#' @export
angelier <- function(x, pch = 1, lwd = 1, lty = 1, col = "black", cex = 1, points = TRUE, bg = NULL, ...) {
  stopifnot(is.Fault(x))

  lines(Fault_plane(x), lwd = lwd, lty = lty, col = col)
  stereo_arrows(Fault_slip(x), sense = x[, "sense"], col = col, ...)
  if (isTRUE(points)) points(Fault_slip(x), pch = pch, col = col, cex = cex, bg = bg)
}



#' Variance visualization
#'
#' Shows the greatcircle of the shortest distance between a set of vectors to a
#' specified vector in a stereoplot.
#' The greatcircles are color-coded by the angular distance.
#'
#' @param x set of vectors. Object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @param y The vector from which the variance should be visualized (only one vector allowed).
#' When `NULL`, then the mean vector of `x` is used (the default).
#' @param .mean character. The type of mean to be used if `y` is `NULL`.
#' One of `"geodesic"` (the default), `"arithmetic"` or `"projected"`.
#' @param ... optional arguments passed to [assign_col()]
#' @param segments logical. Whether the segments should be shown or only the points?
#' @param show.center logical. Whether the center point `y` of the mean should be highlighted in the plot.
#' @inheritParams stereo_point
#'
#' @returns list. `angles` is a vector of the geodesic angles (in degrees)
#' between all vectors in `x` and `y` (or the mean), and `var` is a scalar giving the Fr&eacute;chet variance.
#' @export
#'
#' @seealso [stereo_segment()], [sph_mean()], [geodesic_mean()], [projected_mean()], [geodesic_var()]
#'
#' @examples
#' variance_plot(example_lines)
#' variance_plot(example_planes, example_planes[1, ], segments = FALSE)
variance_plot <- function(x, y = NULL, .mean = c("geodesic", "arithmetic", "projected"), show.center = TRUE, segments = TRUE, upper.hem = FALSE, earea = TRUE, ...) {
  if (is.null(y)) {
    .mean <- match.arg(.mean)
    y <- if (.mean == "geodesic") geodesic_mean(x) else if (.mean == "arithmetic") sph_mean(x) else ot_eigen(x)$vectors[1, ]
  }
  xl <- Line(x)
  yl <- Line(y)

  suppressWarnings(
    ang <- abs(angle(xl, yl))
  )
  ang <- ifelse(ang > 90, 180 - ang, ang)
  cond <- is.nan(ang) | is.na(ang)
  x2 <- xl[!cond, ]

  stereoplot(guides = FALSE, earea = earea)

  if (segments) {
    ang_col <- assign_col(ang[!cond], ...)
    lapply(seq_len(nrow(x2)), function(i) {
      stereo_segment(x2[i, ], yl, col = ang_col[i], upper.hem = upper.hem, earea = earea)
    })
    points(xl, col = "black", pch = 16, cex = .66, upper.hem = upper.hem, earea = earea)
  } else {
    ang_col <- assign_col(ang, ...)
    points(xl, col = ang_col, pch = 16, cex = .66, upper.hem = upper.hem, earea = earea)
  }

  if (isTRUE(show.center)) {
    points(yl, pch = 1, cex = 1, lwd = 2, col = "red", upper.hem = upper.hem, earea = earea)
  }

  var_frechet <- geodesic_var(x, y)

  graphics::title(
    main = "Variance plot",
    sub = paste(
      "Center vector:", round(y[1, 1]), "/", round(y[1, 2]),
      "\nFrechet variance:", round(var_frechet, 2)
    )
  )

  invisible(list(angles = ang, var = var_frechet))
}


#' Plot the bootstrapped confidence ellipse
#'
#' Adds an ellipse marking the bootstrapped confidence interval of the arithmetic mean to an existing plot
#'
#' @param x Spherical object or a list containing the output of an earlier call of [confidence_ellipse()]
#' @param .center logical. Whether the ellipse's center should be plotted?
#' @param col Color of the ellipse and its center
#' @param pch,cex Plotting symbol and size of the ellipse center. Ignored if `.center` is `FALSE`.
#' @param ... graphical parameters passed to [graphics::lines()]
#' @param params list. Parameters passed to [confidence_ellipse()]
#' @inheritParams stereo_smallcircle
#'
#' @returns output of [confidence_ellipse()]
#' @seealso [confidence_ellipse()]
#' @export
#'
#' @examples
#' set.seed(20250411)
#' plot(example_lines, col = "grey")
#' stereo_confidence(example_lines, params = list(n = 100, res = 100), col = "red")
stereo_confidence <- function(x, params = list(), .center = TRUE, col = par("col"), cex = par("cex"), pch = 16, upper.hem = FALSE, earea = TRUE, BALL.radius = 1, ...) {
  if (is.spherical(x)) {
    ce <- do.call(confidence_ellipse, append(list(x = x), params))
  } else if (is.list(x)) {
    ce <- x
  }

  if (.center) {
    e.center <- ce$center
    points(e.center, pch = 16, col = col, cex = cex, upper.hem = upper.hem, earea = earea)
  }

  D <- Line(ce$ellipse)
  Sc <- stereo_coords(D[, 1], D[, 2], upper.hem = upper.hem, earea = earea)

  n <- nrow(D)

  diss <- sqrt((Sc[1:(n - 1), "x"] - Sc[2:(n), "x"])^2 + (Sc[1:(n - 1), "y"] - Sc[2:(n), "y"])^2)
  ww <- which(diss > 0.9 * BALL.radius)
  if (length(ww) > 0) {
    Sc[ww, "x"] <- NA
    Sc[ww, "y"] <- NA
  }
  graphics::lines(Sc[, "x"], Sc[, "y"], col = col, ...)

  invisible(ce)
}
