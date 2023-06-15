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
#' @returns two-column vector with the transformed coordinates
#' @export
stereo_coords <- function(az, inc, upper.hem = FALSE) {
  if (upper.hem) {
    az <- az + 180
  }

  A <- list(az = az, inc = 90 - inc)
  # A$inc <- inc
  # A$az <- az
  B <- fix_inc(A)
  trot <- deg2rad(B$az)
  tinc <- B$inc
  xi <- deg2rad(tinc)
  tq <- sqrt(2) * sin(xi / 2)
  pltx <- tq * sin(trot)
  plty <- tq * cos(trot)

  cbind(x = pltx, y = plty)
}

# stereo2cart <- function(azi, inc, line = TRUE) {
#   deg2rad <- pi / 180
#
#   azir <- azi * deg2rad
#   incr <- inc * deg2rad
#   z <- sin(incr)
#   temp <- cos(incr)
#   x <- cos(azir) * temp
#   y <- sin(azir) * temp
#   return(cbind(x = x, y = y, z = z))
# }


#' Stereographic projection of lines and planes
#'
#' Visualization of lines, planes in a stereographic projection.
#'
#' @param x Object of class `"line"` or `"plane"`
#' @param upper.hem logical. Whether the projection is shown for upper
#' hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
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
stereo_point <- function(x, col = 1, pch = 20, lab = NULL, text.pos = 4, cex = .75, upper.hem = FALSE, ...) {
  stopifnot(is.spherical(x))

  if (is.plane(x) | is.fault(x)) {
    x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  }
  # class(x) <- c("matrix", "array")

  crds <- stereo_coords(
    x[, 1],
    x[, 2],
    upper.hem
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
#' @param ... optional graphical parameters
#' @note `"plane"` objects will be displayed as pole to the plane.
#' @importFrom graphics points text
#' @examples
#' \dontrun{
#' faults <- Fault(
#'   c(0, 90, 180, 270),
#'   c(80, 45, 80, 45),
#'   c(0, 180, 180, 315),
#'   c(80, 5, 80, 36),
#'   c(1, -1, 1, 0)
#' )
#' stereoplot()
#' stereo_fault(faults, col = 1:4)
#' }
stereo_fault <- function(x, hoeppner = FALSE, greatcirles = TRUE, pch = 21, col = 1, lwd = 1, lty = 1, lab = NULL, cex = .75, text.pos = 4, upper.hem = FALSE, ...) {
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
  crds.l <- stereo_coords(
    x[, 3],
    x[, 4],
    upper.hem
  )

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
      points(NULL)
    } else {
      # angle = vangle(Line(x0[i, 3], x0[i, 4]), Line(0, 90)) + 90
      angle <- NULL

      if (x0[i, "sense"] > 0) {
        graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = "\u2191", col = col[i], srt = angle)
      } else if (x0[i, "sense"] < 0) {
        graphics::text(crds.l[i, "x"], crds.l[i, "y"], labels = "\u2191", col = col[i], srt = angle + 180)
      } else {
        graphics::points(crds.l[i, "x"], crds.l[i, "y"], pch = pch[i], col = col[i])
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
#' @param ... optional graphical parameters
#' @importFrom graphics lines
#' @name stereo_cones
#' @examples
#' stereoplot()
#' stereo_point(Line(90, 5), lab = "L")
#' stereo_smallcircle(Line(90, 5), d = 10)
#' stereo_point(Plane(120, 30), lab = "P", col = "red")
#' stereo_greatcircle(Plane(120, 30), col = "red", N = 1000)
NULL

#' @rdname stereo_cones
#' @export
stereo_smallcircle <- function(x, d = 90, col = 1, N = 100, BALL.radius = .25, upper.hem = FALSE, lty = 1, lwd = 1, ...) {
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
    Sc <- stereo_coords(phi2, 90 - theta2, upper.hem)

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
  graphics::lines(c(-1, 1), c(0, 0), col = col)
  graphics::lines(c(0, 0), c(-1, 1), col = col)
}


#' Stereographic projection
#'
#' Initialize the plot for equal-area stereographic projections, i.e. Lambert
#' azimuthal Equal-Area projections (Schmidt).
#'
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
#' @importFrom grDevices gray
#' @importFrom graphics title lines segments
#' @export
#' @examples
#' stereoplot(ticks = 45, title = "title", sub = "subtitle")
stereoplot <- function(guides = TRUE, d = 10, col = grDevices::gray(0.9),
                       lwd = 1, lty = 1, border = "black", title = NULL,
                       sub = NULL, centercross = TRUE, ticks = NULL) {
  plot(c(-1, 1), c(-1, 1),
    type = "n", xlab = NULL, ylab = NULL, asp = 1,
    axes = FALSE, ann = FALSE
  )

  graphics::title(main = title, sub = sub)
  graphics::mtext("N")

  if (guides) stereo_guides(d, col = col, lwd = lwd, lty = lty)

  if (!is.null(ticks)) stereoplot_ticks(angle = ticks, col = border)

  if (centercross) points(0, 0, pch = 3, col = border)
  # graphics::segments(c(-0.02, 0), c(0, -0.02), c(0.02, 0), c(0, 0.02),
  #                    col = "black"
  # )
  stereoplot_frame(col = col, border = border, ndiv = 100)
}

stereoplot_ticks <- function(radius = 1, lenght = 0.02, angle = 10, ...) {
  DR <- radius + lenght
  ang <- pi * seq(0, 360, by = angle) / 180
  segments(
    radius * cos(ang), radius * sin(ang), DR * cos(ang), DR * sin(ang),
    ...
  )
}

stereo_guides <- function(d = 10, col = "black", lwd = 1, lty = 1) {
  lam <- seq(from = 0, to = 180, by = 5) * DEG2RAD()
  lam0 <- pi / 2
  for (j in seq(from = -90 + d, to = 90 - d, by = 10)) {
    phi <- deg2rad(j)
    R <- sqrt(2) / 2
    kp <- sqrt(2 / (1 + cos(phi) * cos(lam - lam0)))
    x <- R * kp * cos(phi) * sin(lam - lam0)
    y <- R * kp * sin(phi)
    lines(x, y, col = col, lwd = lwd, lty = lty)
  }
  phi <- seq(from = -90, to = 90, by = 5) * DEG2RAD()
  for (j in seq(from = d, to = 180 - d, by = d)) {
    lam <- deg2rad(j)
    R <- sqrt(2) / 2
    kp <- sqrt(2 / (1 + cos(phi) * cos(lam - lam0)))
    x <- R * kp * cos(phi) * sin(lam - lam0)
    y <- R * kp * sin(phi)
    lines(x, y, col = col, lwd = lwd, lty = lty)
  }
}
