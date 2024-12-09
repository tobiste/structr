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
  # nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
  nz <- n[, 3]
  # cbind(
  azimuth <- atan2d(n[, 2], n[, 1])
  plunge <- asind(nz)
  # )

  res <- mapply(correct_inc, azi = azimuth, inc = plunge) |> t()
  rownames(res) <- rownames(x)
  colnames(res) <- c("azimuth", "plunge")
  res
}

vec2fol0 <- function(x, y, z) {
  n <- vnorm(cbind(x, y, z)) # normalized vector
  # nz <- sapply(n[, 3], function(x) ifelse(x < 0, -x, x))
  nz <- n[, 3]

  dip_direction <- (atan2d(n[, 2], n[, 1]) + 180)
  dip <- 90 - asind(nz)

  res <- mapply(correct_inc, azi = dip_direction, inc = dip) |> t()
  rownames(res) <- rownames(x)
  colnames(res) <- c("dip_direction", "dip")
  res
}

correct_inc <- function(azi, inc) {
  if (inc > 90) {
    inc <- 180 - 90
    azi <- (azi + 180) %% 360
  } else if (inc < 0) {
    inc <- abs(inc)
    azi <- (azi + 180) %% 360
  }
  c(azi, inc)
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
  if (is.plane(x)) {
    fol2vec0(x[, 1], x[, 2])
  } else {
    lin2vec0(x[, 1], x[, 2])
  }
}

#' @rdname coordinates
#' @export
to_spherical <- function(x, class = c("line", "plane")) {
  class <- match.arg(class)
  if (class == "line") {
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
  structure(
    vec2lin0(x[, 1], x[, 2], x[, 3]),
    class = "line"
  )
}

#' @rdname coordinates
#' @export
vec2plane <- function(x) {
  x <- vec2mat(x)
  structure(vec2fol0(x[, 1], x[, 2], x[, 3]),
    class = "plane"
  )
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
#' `1`or `-1` for normal or thrust offset, respectively.
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
Pair <- function(dip_direction, dip, azimuth, plunge) {
  p <- Fault(dip_direction, dip, azimuth, plunge)
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

# sph2stereo <- function(...){
#   stereo_coords(...)
#   }


#' Converts strike into dip direction using right-hand rule
#' @param strike strike direction in degrees
#' @param dipdirection strike direction in degrees
#'
#' @name rhr
NULL

#' @rdname rhr
#' @export
rhr2dd <- function(strike) {
  (strike + 90) %% 360
}

#' @rdname rhr
#' @export
dd2rhr <- function(dipdirection) {
  (dipdirection - 90) %% 360
}

parse_strike_dip <- function(strike, dip) {
  strike <- parse_azimuth(strike)
  dd <- split_trailing_letters(dip)$measurement
}



#' Parse measurement and direction strings
#'
#' @param x character or number
#'
#' @return list
#' @export
#'
#' @examples
#' test <- c("45NW", "4SE")
#' split_trailing_letters(test)
split_trailing_letters <- function(x) {
  result <- sapply(x, function(x) {
    if (grepl("[NESWnesw]", x)) {
      # Separate numeric and alphabetic parts
      num_part <- as.numeric(gsub("[^0-9.-]", "", x))
      alpha_part <- gsub("[0-9.-]", "", x)
    } else {
      num_part <- as.numeric(x)
      alpha_part <- NA
    }
    list(num = num_part, alpha = alpha_part)
  }, simplify = FALSE)

  # Extract numeric and alpha parts from the list
  number_vector <- sapply(result, function(x) x$num)
  alpha_vector <- sapply(result, function(x) x$alpha)

  # Return as a list
  list(measurement = number_vector, direction = alpha_vector)
}

parse_azimuth <- function(azimuth) {
  sapply(azimuth, function(x) {
    if (is.numeric(x)) {
      parse_quadrant_measurement(x)
    } else {
      stop(
        paste(
          "Ambiguos azimuth:",
          x
        )
      )
    }
  }, simplify = TRUE)
}

rotation_direction <- function(first, second) {
  first_rad <- first * pi / 180
  second_rad <- second * pi / 180
  t(c(cos(first_rad), sin(first_rad))) %*%
    c(cos(second_rad), sin(second_rad))
}

quadrantletter_to_azimuth <- function(x) {
  letters <- trimws(x) |>
    strsplit("") |>
    unlist()
  azimuth <- c("N" = 0, "S" = 180, "E" = 90, "W" = 270)
  tectonicr::circular_mean(azimuth[letters], axial = FALSE)
}


#' Title
#'
#' @param x character.
#'
#' @return numeric
#' @export
#'
#' @examples
#' parse_quadrant_measurement(c("E30N", "W10S"))
parse_quadrant_measurement <- function(x) {
  sapply(x, function(x) {
    x <- trimws(x)

    first_dir <- quadrantletter_to_azimuth(toupper(x[1]))
    sec_dir <- quadrantletter_to_azimuth(toupper(x[length(x)]))

    x2 <- trimws(x) |>
      strsplit("") |>
      unlist()

    angle <- x2[3:length(x2) - 1] |>
      paste(collapse = "") |>
      as.numeric()

    direc <- rotation_direction(first_dir, sec_dir)
    azi <- first_dir + direc * angle

    # Catch ambiguous measurements such as N10S and raise an error
    stopifnot(abs(direc) >= 0.9)
    azi %% 360
  }, simplify = TRUE)
}
