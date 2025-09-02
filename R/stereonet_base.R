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
#' @returns two-column vector with the transformed coordinates
#' @export
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
#' @param x object of class `"Vec3"`, `"Line"`, `"Plane"`. `"Pair"`, or `"Fault"`
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

#' Stereographic projection of faults
#'
#' Visualization of faults (planes and lines) in a stereographic projection.
#'
#' @param x Object of class `"Fault"`
#' @param hoeppner logical. `TRUE` for Hoeppner plot
#' @param greatcirles logical. Whether greatcircles are displayed (`TRUE`, the default) or poles to planes (`FALSE`)
#' @param pch,col,lwd,lty plotting parameters for planes and lines
#' @param lab character. text labels
#' @param cex character expansion of labels
#' @param text.pos position for labels
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param earea logical `TRUE` for Lambert equal-area projection (also "Schmidt net"; the default), or
#' `FALSE` for meridional stereographic projection (also "Wulff net" or "Stereonet").
#' @param ... optional graphical parameters
#' @note `"plane"` objects will be displayed as pole to the plane.
#' @importFrom graphics points text
#'
#' @export
#'
#' @examples
#' faults <- Fault(
#'   c(0, 90, 180, 270),
#'   c(80, 45, 80, 45),
#'   c(0, 170, 180, 315),
#'   c(80, 10, 80, 36),
#'   c(1, -1, 0, 1)
#' )
#' stereoplot()
#' stereo_fault(faults, col = 1:4)
#' # stereo_fault(faults, col =1:4, hoeppner = TRUE)
#' legend("bottomright", c("normal", "thrust", "unknown", "normal"), fill = 1:4)
stereo_fault <- function(x, hoeppner = FALSE, greatcirles = TRUE, pch = 16, col = 1, lwd = 1, lty = 1, lab = NULL, cex = 1, text.pos = 4, upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.Fault(x))
  x0 <- x

  if (length(col) == 1) col <- rep(col, nrow(x))
  if (length(pch) == 1) pch <- rep(pch, nrow(x))
  if (length(cex) == 1) cex <- rep(cex, nrow(x))
  if (length(lwd) == 1) lwd <- rep(lwd, nrow(x))
  if (length(lty) == 1) lty <- rep(lty, nrow(x))

  x[, 1] <- 180 + x[, 1]
  x[, 2] <- 90 - x[, 2]

  crds.p <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem
  )
  if (hoeppner) {
    crds.l <- stereo_coords(
      x[, 3] + 180,
      90 - x[, 4],
      upper.hem, earea
    )
  } else {
    crds.l <- stereo_coords(
      x[, 3],
      x[, 4],
      upper.hem, earea
    )
  }

  for (i in 1:nrow(crds.p)) {
    if (greatcirles) {
      stereo_greatcircle(Plane(x0[i, 1], x0[i, 2]), col = col[i], lwd = lwd[i], lty = lty[i])
    } else {
      graphics::points(crds.p[i, "x"], crds.p[i, "y"], pch = pch[i], col = col[i])
      if (!is.null(lab)) {
        graphics::text(crds.p[, "x"], crds.p[, "y"], labels = lab, pos = text.pos, col = col[i], cex = cex[i])
      }
    }
    if (hoeppner) {
      ang <- (x0[i, 3] * sign(x0[i, "sense"]))
      graphics::points(crds.l[i, "x"], crds.l[i, "y"], pch = pch[i], col = col[i], cex = cex[i] / 1.25)
      if (x0[i, "sense"] != 0) {
        graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = "\u2191", col = col[i], srt = ang, cex = cex[i] * 1.5)
      }
    } else {
      ang <- x0[i, 3] %% 360

      symb <- if (x0[i, 5] < 0) {
        "\u2193"
      } else {
        "\u2191"
      }

      if (x0[i, "sense"] != 0) {
        graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = symb, col = col[i], srt = -ang, cex = cex[i] * 1.5)
      } else {
        graphics::points(crds.l[i, "x"], crds.l[i, "y"], pch = pch[i], col = col[i], cex = cex[i])
      }
    }
  }
}


#' Stereographic projection of cones
#'
#' Visualization of smallcircles and greatcircles in a stereographic projection.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"`, or `"Fault"`.
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


#' Stereoplot frame
#' 
#' Adds a (primitive) circle with given radius to an existing plot
#' 
#' @param ndiv integer. Resolution of circle's line
#' @inheritParams stereoplot
#' @param ... optional arguments passed to [graphics::lines()]
#'
#' @export
#' @seealso [stereo_plot()], [stereoplot_ticks()], [stereoplot_guides()]
#' 
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_frame(col = 'red', lwd = 3)
stereoplot_frame <- function(ndiv = 144, radius = 1, ...) {
  phi <- seq(0, 2 * pi, by = 2 * pi / ndiv)
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
stereoplot <- function(earea = TRUE, guides = TRUE, d = 10, col = grDevices::gray(0.9),
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

  stereoplot_frame(col = border.col, ndiv = 100, radius = radius)
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
#' @seealso [stereo_plot()], [stereoplot_frame()], [stereoplot_guides()]
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_frame()
#' stereoplot_ticks(length = 0.05, angle = 45, col = 'blue', lwd = 2, labels = TRUE)
stereoplot_ticks <- function(length = 0.02, angle = 10, labels = FALSE, ladj = 2 * length, radius = 1, rotation = 0, ...) {
  DR <- radius + length
  ang_deg <- (seq(0, 360, by = angle) + rotation) %% 360
  ang <- pi * ang_deg / 180
  graphics::segments(
    radius * cos(ang), radius * sin(ang), DR * cos(ang), DR * sin(ang),
    ...
  )
  
  do.label <- is.logical(labels)
  if(do.label & isTRUE(labels)){
    labels <- seq(0, 360, by = angle)
    labels[length(labels)] <- NA
    labels <- as.character((-labels+90) %% 360)
    } else labels <- NULL
  if(!is.null(labels)){
    stopifnot(length(labels)==length(ang))
    DR.labs <- DR + ladj
    graphics::text(
      DR.labs * cos(ang), DR.labs * sin(ang), labels = labels, ...
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

  lam_vals <- deg2rad(seq(d, 180 - d, d) )
  
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
#' @seealso [stereo_plot()], [stereoplot_frame()], [stereoplot_ticks()]
#' @export
#' @examples
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_guides(d = 5, earea = FALSE, col = 'green', rotation = 20)
#' 
#' plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
#' stereoplot_guides(d = 15, earea = TRUE, col = 'orange', rotation = 90)
stereoplot_guides <- function(d = 10, earea = TRUE, radius = 1, ...) {
  if (earea) {
    stereo_guides_schmidt(d = d, r = radius, ...)
  } else {
    stereo_guides_wulff(d = d, r = radius, ...)
  }
}


#' Plot spherical objects
#'
#' @param x objects of class `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"`, or `"Fault"`.
#' @inheritParams stereoplot
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param grid.params list.
#' @param ... parameters passed to [stereo_point()], [stereo_smallcircle()], [stereo_greatcircle()], or [stereo_fault()]
#'
#' @exportS3Method graphics::plot 
#'
#' @examples
#' plot(Line(c(90, 80), c(10, 75)), lab = c("L1", "L2"))
#' plot(Plane(120, 30), col = "red")
plot.spherical <- function(x, upper.hem = FALSE, earea = TRUE, grid.params = list(), ...) {
  do.call(stereoplot, append(grid.params, earea))

  if (is.Line(x) | is.Vec3(x)) stereo_point(x, upper.hem = upper.hem, earea = earea, ...)
  if (is.Plane(x)) stereo_greatcircle(x, upper.hem = upper.hem, earea = earea, ...)
  if (is.Fault(x)) stereo_fault(x, upper.hem = upper.hem, earea = earea, ...)
}

# #' @export
# #' @keywords internal
# plot <- function(x, ...) UseMethod("plot", x, ...)

#' Add Points to a Plot
#'
#' @param ... arguments passed to [graphics::points()]
#' @inheritParams plot.spherical
#' @inheritParams graphics::text
#' @importFrom graphics points
#'
#' @exportS3Method graphics::points
#'  
#' @examples
#' stereoplot()
#' points(rvmf(n = 100))
#'
#' points(Plane(120, 30), col = "red")
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

# #' @export
# #' @keywords internal
# points <- function(x, ...) UseMethod("points", x, ...)


#' Add Lines to a Plot
#'
#' @inheritParams plot.spherical
#' @inheritParams graphics::text
#' @param ... arguments passed to [graphics::lines()]
#' @importFrom graphics lines
#' @exportS3Method graphics::lines
#'
#' @examples
#' stereoplot()
#' lines(rvmf(n = 5), d = runif(5, 0, 90), col = 1:5)
lines.spherical <- function(x, ...) stereo_smallcircle(x, ...)

# #' @export
# #' @keywords internal
# lines <- function(x, ...) UseMethod("lines", x, ...)


#' Add Points to a Plot
#'
#' @param ... arguments passed to [graphics::text()]
#' @inheritParams plot.spherical
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

# #' @export
# #' @keywords internal
# text <- function(x, ...) UseMethod("text", x, ...)
