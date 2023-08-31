fol2vec0 <- function(azi, inc) {
  azi <- deg2rad(azi)
  inc <- deg2rad(inc)
  cbind(
    x = -cos(azi) * sin(inc),
    y = -sin(azi) * sin(inc),
    z = cos(inc)
  )
}

lin2vec0 <- function(azi, inc) {
  azi <- deg2rad(azi)
  inc <- deg2rad(inc)
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
    azimuth = atan2d(n[, 2], n[, 1]) %% 360,
    plunge = asind(nz)
  )
}

vec2fol0 <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
  cbind(
    dip_direction = (atan2d(n[, 2], n[, 1]) + 180) %% 360,
    dip = 90 - asind(nz)
  )
}


#' Cartesian and spherical coordinates
#'
#' Converts between Cartesian and spherical coordinates
#'
#' @param l,p two column array, or two element vector. Spherical coordinates
#' (in degrees) of lines (`"azimuth"`, `"plunge"`) or plane (`"dip_direction"`,
#' `"dip"`).
#' @param x 3three-column array, or three-element vector. Cartesian coordinates
#' @param class character. Either `"line"` or `"plane"`
#' @name coordinates
#' @examples
#' cbind(c(90, 89), c(45, 46)) |>
#'   line2vec() |>
#'   vec2line()
NULL

#' @rdname coordinates
#' @export
to_vec <- function(x) {
  # x <- vec2mat(x)
  if(is.plane(x)){
    fol2vec0(x[, 1], x[, 2])
  } else {
    lin2vec0(x[, 1], x[, 2])
  }
}

#' @rdname coordinates
#' @export
to_spherical <- function(x, class = c("line", "plane")) {
  class <- match.arg(class)
  if (is.line(x)) {
    vec2line(x)
  } else {
    vec2plane(x)
  }
}

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
vec2line <- function(x) {
  x <- vec2mat(x)
  l <- vec2lin0(x[, 1], x[, 2], x[, 3])
  class(l) <- "line"
  l
}

#' @rdname coordinates
#' @export
vec2plane <- function(x) {
  x <- vec2mat(x)
  p <- vec2fol0(x[, 1], x[, 2], x[, 3])
  class(p) <- "plane"
  p
}


geo2vec0 <- function(lat, lon) {
  cbind(
    x = cosd(lat) * cosd(lon),
    y = cosd(lat) * sind(lon),
    z = sind(lat)
  )
}

geo2vec <- function(g) {
  x <- vec2mat(g)
  geo2vec0(x[, 1], x[, 2])
}


#' Structure classes
#'
#' @description
#' `Line`, `Plane`, and `Fault` create a `"line"`, `"plane"`, `"pair"`, and `"fault"`
#' class object, respectively, from the given set of values.
#'
#' `as.line`, `as.plane`, `as.pair`, and `as.fault` attempt to turn its argument into a
#' `"line"`, `"plane"`, and `"fault"` class object, respectively.
#'
#' `is.line`, `is.plane`, `is.pair`, and `is.fault` test if its argument is a `"line"`,
#' `"plane"`, and `"fault"` class object, respectively.
#'
#' @param l,p,f numeric vector or array containing the spherical coordinates
#' (1st element/column is azimuth, 2nd element/column is inclination, both in
#' degrees), or object of class `"line"`, `"plane"`, `"pair"`, or `"fault"`
#' @param azimuth,plunge numeric vectors. Azimuth and plunge of a line (in
#' degrees)
#' @param dip_direction,dip numeric vectors. Dip direction and dip of a plane
#' (in degrees)
#' @param sense (optional) integer. Sense of the line on a fault plane. Either
#' `1`or `-1` for normal/dextral or thrust/sinistral offset, respectively.
#' @details
#' `is.line`, `is.plane`, `is.pair`, and `is.fault` return `TRUE` if `l`, `p`, and `f`
#' are an object of class `"line"`, `"plane"`, `"pair"` or `"fault"`, respectively, and
#' `FALSE` otherwise.
#'
#' `is.spherical` returns `TRUE` if the argument's class is one of `"line"`,
#' `"plane"`, `"pair"`, or `"fault"` and `FALSE` otherwise
#'
#' `as.line`, `as.plane`, `as.pair`, and `as.fault` are is generic functions. If the
#' argument is a `"line"` or `"plane"` class, it will be converted.
#'
#' @name classes
#' @examples
#' x <- Line(120, 50) # create line
#' is.line(x) # test if line
#' as.plane(x) # convert to plane
NULL


#' @rdname classes
#' @export
Line <- function(azimuth, plunge) {
  stopifnot(length(azimuth) == length(plunge))
  cbind(azimuth, plunge) |> as.line()
}

#' @rdname classes
#' @export
Plane <- function(dip_direction, dip) {
  stopifnot(length(dip_direction) == length(dip))
  cbind(dip_direction, dip) |> as.plane()
}

#' @rdname classes
#' @export
Fault <- function(dip_direction, dip, azimuth, plunge, sense = NULL) {
  stopifnot(
    length(dip_direction) == length(dip),
    length(dip_direction) == length(azimuth),
    length(dip_direction) == length(plunge)
  )
  if (is.null(sense)) {
    sense <- rep(0, length(dip_direction))
  } else {
    stopifnot(
      length(dip_direction) == length(sense)
    )
  }
  cbind(dip_direction, dip, azimuth, plunge, sense) |> as.fault()
}

#' @rdname classes
#' @export
Pair <- function(dip_direction, dip, azimuth, plunge){
  p = Fault(dip_direction, dip, azimuth, plunge) 
  as.pair(p[, -5])
}

#' @rdname classes
#' @export
as.line <- function(l) {
  if (is.plane(l)) {
    x <- l |>
      plane2vec() |>
      vec2line()
  } else {
    x <- vec2mat(l)
    colnames(x) <- c("azimuth", "plunge")
  }
  class(x) <- "line"
  rownames(x) <- rownames(l)
  x
}

#' @rdname classes
#' @export
as.plane <- function(p) {
  if (is.line(p)) {
    x <- p |>
      line2vec() |>
      vec2plane()
  } else {
    x <- vec2mat(p)
    colnames(x) <- c("dip_direction", "dip")
  }
  class(x) <- "plane"
  rownames(x) <- rownames(p)
  x
}

#' @rdname classes
#' @export
as.fault <- function(f) {
  x <- vec2mat(f)
  class(x) <- c("fault", "pair")
  colnames(x) <- c("dip_direction", "dip", "azimuth", "plunge", "sense")
  rownames(x) <- rownames(f)
  x
}

#' @rdname classes
#' @export
as.pair <- function(f) {
  x <- vec2mat(f)
  class(x) <- "pair"
  colnames(x) <- c("dip_direction", "dip", "azimuth", "plunge")
  rownames(x) <- rownames(f)
  x
}

#' @rdname classes
#' @export
is.line <- function(l) {
  inherits(l, "line")
}

#' @rdname classes
#' @export
is.plane <- function(p) {
  inherits(p, "plane")
}

#' @rdname classes
#' @export
is.fault <- function(f) {
  inherits(f, "fault")
}

#' @rdname classes
#' @export
is.pair <- function(f) {
  inherits(f, "pair")
}

#' @rdname classes
#' @export
is.spherical <- function(l) {
  inherits(l, "plane") | inherits(l, "line") | inherits(l, "pair")
}

#' Spherical coordinate conversion
#'
#' @inheritParams best_cone_ramsay
#' @details acoscartesian coordinates are given as
#' \describe{
#' \item{\code{a}}{deviation from x-axis E-W horizontal (angle in degrees)}
#' \item{\code{b}}{deviation from y-axis N-S horizontal (angle in degrees)}
#' \item{\code{c}}{deviation from z-axis vertical (angle in degrees)}
#' }
#' @source Ramsay, 1967, p. 15-16
#' @name ramsay_coords
#' @importFrom tectonicr deg2rad rad2deg
#' @examples
#' \dontrun{
#' # Stereographic coordinates (angle notation):
#' x <- rbind(
#'   c(58, 78, 35),
#'   c(47, 72, 47),
#'   c(57, 68, 41),
#'   c(56, 64, 45),
#'   c(52, 65, 47),
#'   c(61, 62, 51),
#'   c(55, 60, 49),
#'   c(61, 62, 42),
#'   c(68, 52, 46),
#'   c(42, 58, 65)
#' )
#' acoscartesian_to_cartesian(x)
#' acoscartesian_to_cartesian(x) |> to_spherical()
#' acoscartesian_to_cartesian(x) |> cartesian_to_acoscartesian()
#' }
NULL

#' @rdname ramsay_coords
#' @export
cartesian_to_acosvec <- function(x) {
  # acos(x) |> tectonicr::rad2deg()
  a <- acos(x[, 1])
  b <- acos(x[, 2])
  c <- acos(x[, 3])

  a <- ifelse(x[, 1] < 0, -a, a)
  b <- ifelse(x[, 2] < 0, -b, b)
  c <- ifelse(x[, 3] < 0, -c, c)
  tectonicr::rad2deg(cbind(a, b, c))
}

#' @rdname ramsay_coords
#' @export
acoscartesian_to_cartesian <- function(x) {
  # tectonicr::deg2rad(x) |> cos()
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1])
  cy <- cos(x[, 2])
  cz <- cos(x[, 3])

  cx <- ifelse(x[, 1] < 0, -cx, cx)
  cy <- ifelse(x[, 2] < 0, -cy, cy)
  cz <- ifelse(x[, 3] < 0, -cz, cz)
  cbind(x = cx, y = cy, z = cz)
}
