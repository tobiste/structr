DEG2RAD <- function() {
  pi / 180
}

roty3 <- function(deg) {
  rad1 <- deg * pi / 180
  r <- diag(3)
  r[1, 1] <- cos(rad1)
  r[3, 1] <- sin(rad1)
  r[3, 3] <- r[1, 1]
  r[1, 3] <- -r[3, 1]
  return(r)
}

rotz3 <- function(deg) {
  rad1 <- deg * pi / 180
  r <- diag(3)
  r[1, 1] <- cos(rad1)
  r[1, 2] <- sin(rad1)
  r[2, 2] <- r[1, 1]
  r[2, 1] <- -r[1, 2]
  return(r)
}

FixDip <- function(A) {
  az <- A$az
  dip <- A$dip
  tdip <- DEG2RAD() * (dip %% 360)
  co <- cos(tdip)
  si <- sin(tdip)
  ang <- atan2(si, co) / DEG2RAD()
  quad <- rep(1, length(dip))
  quad[co >= 0 & si >= 0] <- 1
  quad[co < 0 & si >= 0] <- 2
  quad[co < 0 & si < 0] <- 3
  quad[co >= 0 & si < 0] <- 4
  dip[quad == 1] <- ang[quad == 1]
  dip[quad == 2] <- 180 - ang[quad == 2]
  dip[quad == 3] <- 180 + ang[quad == 3]
  dip[quad == 4] <- -ang[quad == 4]
  az[quad == 1] <- az[quad == 1]
  az[quad == 2] <- 180 + az[quad == 2]
  az[quad == 3] <- az[quad == 3]
  az[quad == 4] <- 180 + az[quad == 4]
  A$az <- az %% 360
  A$dip <- dip
  return(A)
}

stereo_coords <- function(az, iang) {
  sph <- cos(DEG2RAD() * iang)
  sph[sph >= 0] <- 0
  sph[sph < 0] <- 1



  A <- list(az = az, dip = iang)
  A$dip <- iang
  A$az <- az
  B <- FixDip(A)
  trot <- DEG2RAD() * B$az
  tdip <- B$dip
  xi <- DEG2RAD() * tdip
  tq <- sqrt(2) * sin(xi / 2)
  pltx <- tq * sin(trot)
  plty <- tq * cos(trot)

  cbind(x = pltx, y = plty)
}


#' Plot lines and planes
#'
#' Visualization of lines, planes in a stereographic projection.
#' @param azi,inc numeric. Angles of azimuth (or dip direction) and inclination (dip) in degrees.
#' @param plane logical. Whether the object is a linear (`FALSE`) of planar feature (`TRUE`).
#' @param upper.hem logical. Whether the projection is shown for upper hemisphere (`TRUE`) or lower hemisphere (`FALSE`, the default).
#' @param col color
#' @param pch plotting character
#' @param lab character. text labels
#' @param text.pos position for labels
#' @param cex character expansion of labels
#' @param ... optional graphical parameters
#' @export
#' @examples
#' stereoplot()
#' stereo_point(90, 80, lab = "L")
#' stereo_point(120, 30, plane = TRUE, lab = "P", col = "red")
stereo_point <- function(azi, inc, plane = FALSE, upper.hem = FALSE, col = 1, pch = 20, lab = NULL, text.pos = 4, cex = .75, ...) {
  crds <- stereo_coords(
    ifelse(upper.hem, azi + 180, azi),
    ifelse(plane, -inc, 90 - inc)
  )


  points(crds[, "x"], crds[, "y"], pch = pch, col = col, cex = cex, ...)
  if (!is.null(lab)) {
    text(crds[, "x"], crds[, "y"], labels = lab, pos = text.pos, col = col)
  }
}




#' Plot cones
#'
#' Visualization of smallcircles and greatcircles in a stereographic projection.
#'
#' @param azi,inc numeric. Angles of azimuth (or dip direction) and inclination (dip) in degrees.
#' @param d numeric. conical angle in degrees.
#' @param col color
#' @param N integer. number of points to calculate
#' @param BALL.radius numeric size of sphere
#' @param ... optional graphical parameters
#' @export
#' @examples
#' stereoplot()
#' stereo_point(90, 80, lab = "L")
#' stereo_smallcircle(90, 80, d = 10)
#' stereo_point(120, 30, plane = TRUE, lab = "P", col = "red")
#' stereo_greatcircle(120, 30, col = "red", N = 1000)
stereo_smallcircle <- function(azi, inc, d = 90, col = 1, N = 100, BALL.radius = .25, ...) {
  az <- azi
  iang <- 90 - inc
  phi <- seq(from = 0, to = 2 * pi, length = N)
  theta <- d * DEG2RAD()
  x <- BALL.radius * sin(theta) * cos(phi)
  y <- BALL.radius * sin(theta) * sin(phi)
  z <- rep(BALL.radius * cos(theta), length(phi))
  D <- cbind(x, y, z)
  ry <- roty3(iang)
  rz <- rotz3(az)
  Rmat <- ry %*% rz
  g <- D %*% Rmat
  r2 <- sqrt(g[, 1]^2 + g[, 2]^2 + g[, 3]^2)
  phi2 <- atan2(g[, 2], g[, 1])
  theta2 <- acos(g[, 3] / r2)
  Sc <- stereo_coords(phi2 / DEG2RAD(), theta2 / DEG2RAD())

  diss <- sqrt((Sc[1:(N - 1), "x"] - Sc[2:(N), "x"])^2 + (Sc[1:(N - 1), "y"] - Sc[2:(N), "y"])^2)
  ww <- which(diss > 0.9 * BALL.radius)
  if (length(ww) > 0) {
    Sc[ww, "x"] <- NA
    Sc[ww, "y"] <- NA
  }

  lines(Sc[, "x"], Sc[, "y"], col = col, ...)
}

stereo_greatcircle <- function(azi, inc, col = 1, lwd = 1, N = 100, ...) {
  stereo_smallcircle(azi, 90 + inc, d = 90, col = col, lwd = lwd, N = N, ...) # add circle
}


#' Circle plot
pcirc <- function(col = "black", border = "black", ndiv = 36) {
  phi <- seq(0, 2 * pi, by = 2 * pi / ndiv)
  x <- cos(phi)
  y <- sin(phi)
  lines(x, y, col = border)
  lines(c(-1, 1), c(0, 0), col = col)
  lines(c(0, 0), c(-1, 1), col = col)
}


#' Initialize the Stereo plot
#'
#' Plot equal-area stereographic projection, i.e. Lambert azimuthal Equal-Area projection (Schmidt).
#' @param add logical, `TRUE` adds to existing plot
#' @param col	color of lines
#' @param border color of outer rim of stereo plot
#' @param lwd	linewidth of lines description
#' @source [RFOC::net()]
#' @export
#' @examples
#' stereoplot()
stereoplot <- function(add = FALSE, col = gray(0.9), border = "black", lwd = 1, main = "N") {
  plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE, ann = FALSE)
  title(main)

  lam <- seq(from = 0, to = 180, by = 5) * DEG2RAD()
  lam0 <- pi / 2
  for (j in seq(from = -80, to = 80, by = 10)) {
    phi <- j * DEG2RAD()
    R <- sqrt(2) / 2
    kp <- sqrt(2 / (1 + cos(phi) * cos(lam - lam0)))
    x <- R * kp * cos(phi) * sin(lam - lam0)
    y <- R * kp * sin(phi)
    lines(x, y, col = col, lwd = lwd)
  }
  phi <- seq(from = -90, to = 90, by = 5) * DEG2RAD()
  for (j in seq(from = 10, to = 170, by = 10)) {
    lam <- j * DEG2RAD()
    R <- sqrt(2) / 2
    kp <- sqrt(2 / (1 + cos(phi) * cos(lam - lam0)))
    x <- R * kp * cos(phi) * sin(lam - lam0)
    y <- R * kp * sin(phi)
    lines(x, y, col = col, lwd = lwd)
  }
  segments(c(-0.02, 0), c(0, -0.02), c(0.02, 0), c(0, 0.02),
    col = "black"
  )
  pcirc(col = col, border = border, ndiv = 72)
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
