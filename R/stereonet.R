roty3 <- function(deg) {
  rad1 <- deg2rad(deg)
  r <- diag(3)
  r[1, 1] <- cos(rad1)
  r[3, 1] <- sin(rad1)
  r[3, 3] <- r[1, 1]
  r[1, 3] <- -r[3, 1]
  return(r)
}

rotz3 <- function(deg) {
  rad1 <- deg2rad(deg)
  r <- diag(3)
  r[1, 1] <- cos(rad1)
  r[1, 2] <- sin(rad1)
  r[2, 2] <- r[1, 1]
  r[2, 1] <- -r[1, 2]
  return(r)
}

fix_inc <- function(A) {
  az <- A$az
  inc <- A$inc
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
  A$az <- az %% 360
  A$inc <- tinc
  return(A)
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

 
  A <- list(az = az, inc = 90 - inc)
  B <- fix_inc(A)
  azi <- deg2rad(B$az)
  inc <- deg2rad(B$inc)

  if (earea) {
    tq <- sqrt(2) * sin(inc / 2)
  } else {
    tq <- tan(inc / 2)
  }
  pltx <- r * tq * sin(azi)
  plty <- r * tq * cos(azi)
  cbind(x = pltx, y = plty)
}

.schmidt_crds <- function(az_rad, inc_rad, r = 1) {
  tq <- r * sqrt(2) * sin(pi / 4 - inc / 2)
  x <- tq * sin(azi)
  y <- tq * cos(azi)
  cbind(x = x, y = y)
}

.wulff_crds <- function(az_rad, inc_rad, r = 1) {
  tq <- r * tan(pi / 4 - inc / 2)
  x <- tq * sin(azi)
  y <- tq * cos(azi)
  cbind(x = x, y = y)
}


#' Stereographic projection of lines and planes
#'
#' Visualization of lines, planes in a stereographic projection.
#'
#' @param x Object of class `"line"` or `"plane"`
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
#' @note `"plane"` and `"fault"` objects will be displayed as pole to the plane (only).
#' @importFrom graphics points text
#' @export
#' @examples
#' stereoplot()
#' stereo_point(Line(azimuth = c(90, 80), plunge = c(10, 75)), lab = c("L1", "L2"))
#' stereo_point(Plane(120, 30), lab = "P", col = "red")
stereo_point <- function(x, col = 1, pch = 20, lab = NULL, text.pos = 4, cex = 1, upper.hem = FALSE, earea = TRUE, ...) {
  stopifnot(is.spherical(x))

  if (is.plane(x) | is.fault(x)) {
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
#' @param x Object of class `"fault"`
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
  stopifnot(is.fault(x))
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
      # ang = ifelse(ang > 180, 360 - ang, ang)
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
#' @param x Object of class `"line"` or `"plane"`
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
NULL

#' @rdname stereo_cones
#' @export
stereo_smallcircle <- function(x, d = 90, col = 1, N = 1000, upper.hem = FALSE, earea = TRUE, lty = 1, lwd = 1, BALL.radius = 1, ...) {
  stopifnot(is.spherical(x))
  if (length(col) == 1) col <- rep(col, nrow(x))
  if (length(lty) == 1) lty <- rep(lty, nrow(x))
  if (length(lwd) == 1) lwd <- rep(lwd, nrow(x))

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


#' Circle plot
#' @noRd
stereoplot_frame <- function(col = "black", border = "black", ndiv = 36) {
  phi <- seq(0, 2 * pi, by = 2 * pi / ndiv)
  x <- cos(phi)
  y <- sin(phi)
  graphics::lines(x, y, col = border)
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
#' @param border color of outer rim of stereo plot
#' @param title,sub character. Title and subtitle of the plot
#' @param centercross logical. Whether a center cross should be added (`TRUE` by default)
#' @param ticks integer. Angle between ticks. if `NULL` (the default), no ticks are drawn.
#' @source Adapted from the `RFOC` package
#' @importFrom graphics points
#' @export
#' @examples
#' stereoplot(ticks = 45, title = "title", sub = "subtitle")
stereoplot <- function(earea = TRUE, guides = TRUE, d = 10, col = grDevices::gray(0.9),
                       lwd = 1, lty = 1, border = "black", title = NULL,
                       sub = NULL, centercross = TRUE, ticks = NULL) {
  plot(c(-1, 1), c(-1, 1),
    type = "n", xlab = NULL, ylab = NULL, asp = 1,
    axes = FALSE, ann = FALSE
  )

  graphics::title(main = title, sub = sub)
  graphics::mtext("N")

  if (guides) stereo_guides(d=d, earea = earea, col = col, lwd = lwd, lty = lty)

  if (!is.null(ticks)) stereoplot_ticks(angle = ticks, col = border)

  if (centercross) points(0, 0, pch = 3, col = border)

  stereoplot_frame(col = col, border = border, ndiv = 100)
}

#' @importFrom graphics segments
stereoplot_ticks <- function(radius = 1, length = 0.02, angle = 10, ...) {
  DR <- radius + length
  ang <- pi * seq(0, 360, by = angle) / 180
  segments(
    radius * cos(ang), radius * sin(ang), DR * cos(ang), DR * sin(ang),
    ...
  )
}

stereo_guides_schmidt <- function(d = 10, n = 512, ...) {
  lam_seq <- seq(0, 180, length.out = n) * pi / 180
  lam0 <- pi / 2
  R <- sqrt(2) / 2
  
  # Precompute sin and cos of lam - lam0
  cos_lam <- cos(lam_seq - lam0)
  sin_lam <- sin(lam_seq - lam0)
  
  # Latitude lines (constant phi)
  phi_vals <- seq(-90 + d, 90 - d, 10) * pi / 180
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
  phi_seq <- seq(-90, 90, 5) * DEG2RAD()
  cos_phi_seq <- cos(phi_seq)
  sin_phi_seq <- sin(phi_seq)
  
  lam_vals <- seq(d, 180 - d, d) * DEG2RAD()
  cos_lam_vals <- cos(lam_vals - lam0)
  sin_lam_vals <- sin(lam_vals - lam0)
  
  for (j in seq_along(lam_vals)) {
    lam <- lam_vals[j]
    kp <- sqrt(2 / (1 + cos_phi_seq * cos_lam_vals[j]))
    x <- R * kp * cos_phi_seq * sin_lam_vals[j]
    y <- R * kp * sin_phi_seq
    graphics::lines(x, y, ...)
  }
}

stereo_guides_wulff <- function(d = 9, r = 1, n = 512, ...) {
  if (n %% 2 != 0) n <- n + 1
  beta0 <- seq(0, 360, length.out = n) * DEG2RAD()
  beta <- c(beta0[(n / 2):n], beta0[1:(n / 2 - 1)])
  
  cos_beta <- cos(beta)
  sin_beta <- sin(beta)
  
  # great circles
  phi <- seq(0, 180, by = d) * DEG2RAD()
  rcos_phi <- r / cos(phi)
  rtan_phi <- r * tan(phi)
  
  
  for (i in seq_along(phi)) {
    xg <- -rtan_phi[i] + (rcos_phi[i]) * cos_beta
    yg <- rcos_phi[i] * sin_beta
    condg <- sqrt(xg^2 + yg^2) <= r
    graphics::lines(xg[condg], yg[condg], ...)
  }
  
  
  # small circles
  gamma <- seq(0, 180, by = d) * DEG2RAD()
  rtan_gamma <- r * tan(gamma)
  rcos_gamma <- r / cos(gamma)
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


#' @importFrom graphics lines
stereo_guides <- function(d = 10, earea = TRUE, ...) {
  if(earea){
    stereo_guides_schmidt(d = d, ...)
  } else {
    stereo_guides_wulff(d = d, ...)
  }
}

