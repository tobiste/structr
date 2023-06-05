fol2vec0 <- function(azi, inc) {
  azi <- azi * DEG2RAD()
  inc <- inc * DEG2RAD()
  cbind(
    x = -cos(azi) * sin(inc),
    y = -sin(azi) * sin(inc),
    z = cos(inc)
  )
}

lin2vec0 <- function(azi, inc) {
  azi <- azi * DEG2RAD()
  inc <- inc * DEG2RAD()
  cbind(
    x = cos(azi) * cos(inc),
    y = sin(azi) * cos(inc),
    z = sin(inc)
  )
}

vec2lin0 <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
  cbind(
    azi = (atan2(n[, 2], n[, 1]) / DEG2RAD()) %% 360,
    inc = asin(nz) / DEG2RAD()
  )
}

vec2fol0 <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
  cbind(
    azi = ((atan2(n[, 2], n[, 1]) / DEG2RAD()) + 180) %% 360,
    inc = 90 - (asin(nz) / DEG2RAD())
  )
}


#' Coordinate conversion
#'
#' @param l,p two column matrix, or two element vector. Spherical coordinates
#' (in degrees) of lines (azimuth, plunge) or plane (dir direction, dip).
#' @param v three column matrix, or three-element vector. Cartesian coordinates
#' @name coordinates
#' @examples
#' cbind(c(90, 89), c(45, 46)) |>
#'   lin2vec() |>
#'   vec2lin()
NULL

#' @rdname coordinates
#' @export
line2vec <- function(l) {
  x <- vec2mat(l)
  lin2vec0(x[, 1], x[, 2])
}

#' @rdname coordinates
#' @export
plane2vec <- function(p) {
  x <- vec2mat(p)
  fol2vec0(x[, 1], x[, 2])
}

#' @rdname coordinates
#' @export
vec2line <- function(v) {
  x <- vec2mat(v)
  vec2lin0(x[, 1], x[, 2], x[, 3])
}

#' @rdname coordinates
#' @export
vec2plane <- function(v) {
  x <- vec2mat(v)
  vec2fol0(x[, 1], x[, 2], x[, 3])
}
