deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

sind <- function(x) sin(deg2rad(x))
cosd <- function(x) cos(deg2rad(x))
tand <- function(x) tan(deg2rad(x))
asind <- function(x) rad2deg(asin(x))
acosd <- function(x) rad2deg(acos(x))
atand <- function(x) rad2deg(atan(x))
atan2d <- function(y, x) rad2deg(atan2(y, x)) %% 360




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

#' @keywords internal
#' @noRd
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

#' @keywords internal
#' @noRd
parse_azimuth <- function(azimuth) {
  sapply(azimuth, function(x) {
    if (is.numeric(x)) {
      parse_quadrant_measurement(x)
    } else {
      stop(
        paste(
          "Ambiguous azimuth:",
          x
        )
      )
    }
  }, simplify = TRUE)
}

#' @keywords internal
#' @noRd
rotation_direction <- function(first, second) {
  first_rad <- first * pi / 180
  second_rad <- second * pi / 180
  t(c(cos(first_rad), sin(first_rad))) %*%
    c(cos(second_rad), sin(second_rad))
}

#' @keywords internal
#' @noRd
quadrantletter_to_azimuth <- function(x) {
  letters <- trimws(x) |>
    strsplit("") |>
    unlist()
  azimuth <- c("N" = 0, "S" = 180, "E" = 90, "W" = 270)
  tectonicr::circular_mean(azimuth[letters], axial = FALSE)
}

#' Converts azimuth angles into Cardinal directions
#'
#' @param x angles in degree.
#' @param n_directions either 8 for 8-point (N, NE, E, …) or 6 for 16-point (N, NNE, NE, …) cardinal version.
#'
#' @returns character vector
#' @export
#'
#' @examples
#' azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350)) # 8-point compass
#' azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 16) # 16-point compass
azimuth_to_cardinal <- function(x, n_directions = 8) {
  # Normalize to 0–360
  azimuth <- x %% 360

  if (n_directions == 8) {
    dirs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  } else if (n_directions == 16) {
    dirs <- c(
      "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
      "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
    )
  } else {
    stop("n_directions must be 8 or 16")
  }

  step <- 360 / n_directions
  # Round to nearest sector center
  idx <- floor((azimuth + step / 2) / step) + 1
  idx[idx > n_directions] <- 1

  dirs[idx]
}


#' Quadrant measurement expressions to angles
#'
#' Interprets quadrant measurement expressions, such as "E30N" or "W10S" as azimuth angles
#'
#' @param x character.
#'
#' @return numeric
#' @noRd
#' @keywords internal
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


.scale <- function(x, from = range(x), to) {
  original_min <- from[1]
  original_max <- from[2]

  target_min <- to[1]
  target_max <- to[2]

  target_min + (x - original_min) * (target_max - target_min) / (original_max - original_min)
}


# Color assignment helper functions --------------------------------------------

.normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' Helper functions for color assignment and legends
#'
#' @param x vector to colorize
#' @param breaks numeric.
#' @param n integer.
#' @param title character. Legend title
#' @param pal color function; Default is [viridis::viridis()]
#' @param fill color vector
#' @param labels character.vector. Names of discrete colors.
#' Can be ignored when `cols` is a named vector.
#' @param position Legend position
#' @param ... arguments passed to color function
#'
#' @return .
#' @import viridis
#' @importFrom grDevices colorRamp
#' @name colorize
#'
#' @examples
#' set.seed(1234)
#'
#' # example for discrete colors
#' x <- rvmf(5, mu = Line(120, 50), k = 5)
#' key <- letters[round(runif(5, 1, 26))]
#' plot(x, col = assign_col_d(key), grid.params = list(guides = FALSE))
#' legend_d(assign_col_d(key))
#'
#' # example for continuous colors:
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' plot(x, col = assign_col(runif(100)), grid.params = list(guides = FALSE))
#' legend_c(seq(0, 1, .1), title = "test")
NULL

#' @rdname colorize
#' @export
assign_col_d <- function(x, pal = viridis::viridis, ...) {
  groups <- unique(x)
  n <- length(groups)
  cols <- do.call(pal, args = list(n = n, ...))
  named_cols <- stats::setNames(cols, groups)
  named_cols[x]
}

#' @rdname colorize
#' @export
assign_col <- function(x, n = length(x), pal = viridis::viridis, ...) {
  normalized_data <- .normalize(x)
  colors <- do.call(pal, args = list(n = n, ...))
  colors[as.numeric(cut(normalized_data, breaks = n))]
}

#' @rdname colorize
#' @export
assign_col_binned <- function(x, breaks, pal = viridis::viridis, ...) {
  breaks <- pretty(x, n = breaks)
  n2 <- length(breaks) - 1
  cols <- do.call(pal, args = list(n2, ...)) # [order]
  named_cols <- cut(x, breaks = breaks, labels = cols, include.lowest = TRUE) |>
    as.character()
  names(named_cols) <- cut(x, breaks = breaks, include.lowest = TRUE)
  named_cols
}

#' @importFrom grDevices colorRamp
color_func <- function(x, pal = viridis::viridis, ...) {
  color_func0 <- colorRamp(do.call(pal, args = list(n = 10000, ...)))

  grDevices::rgb(color_func0(x, ...) / 255)
}

#' @rdname colorize
#' @export
legend_c <- function(breaks, title = NULL, pal = viridis::viridis, ...) {
  label_pos <- .normalize(breaks)
  legend_image <- grDevices::as.raster(
    matrix(
      color_func(seq(0, 1, .001), pal = pal, ...),
      ncol = 1
    )
  )

  graphics::par(new = TRUE)
  graphics::layout(matrix(1, 1))

  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
  graphics::text(x = 0.9 + (0.1 * 1.1), y = label_pos, labels = breaks, adj = 0)
  graphics::text(x = .925, y = .5, labels = title, adj = 0.5, srt = 90, font = 2)
  graphics::rasterImage(legend_image, .95, 0, 1, 1)
}

#' @rdname colorize
#' @export
legend_d <- function(fill, labels = names(fill), position = "topright", ...) {
  graphics::legend(position,
    legend = labels,
    fill = fill,
    ...
  )
}

#' Add an Ellipse to existing plot
#'
#' @param x,y the x and y co-ordinates for the center(s) of the ellipse(s).
#' @param radius.x  scalar or a vector giving the semi-major axis of the ellipse.
#' @param radius.y a scalar or a vector giving the semi-minor axis of the ellipse.
#' @param rot angle of rotation in radians.
#' @param nv number of vertices to draw the ellipses.
#' @param border color for borders. The default is par("fg"). Use border = NA to omit borders.
#' @param col color(s) to fill or shade the annulus sector with. The default NA (or also NULL) means do not fill (say draw transparent).
#' @param lty line type for borders and shading; defaults to "solid".
#' @param lwd line width for borders and shading.
#' @param plot logical. If TRUE the structure will be plotted. If FALSE only the points are calculated and returned. Use this if you want to combine several geometric structures to a single polygon.
#'
#' @returns The function invisibly returns a list of the calculated coordinates for all shapes.
#' @importFrom graphics par polygon
#' @export
#' @examples
#' plot(c(0, 1), c(0, 1), type = "n")
#' ellipse(.5, .5, radius.x = 0.5, radius.y = .25, col = "darkgreen", border = "red")
ellipse <- function(
    x = 0, y = x,
    radius.x = 1, radius.y = radius.x,
    rot = 0, nv = 512,
    border = par("fg"), col = par("bg"),
    lty = par("lty"), lwd = par("lwd"),
    plot = TRUE) {
  # normalize inputs to same length
  lgp <- list(
    x = x, y = y, radius.x = radius.x,
    radius.y = radius.y, rot = rot, nv = nv
  )
  maxdim <- max(vapply(lgp, length, 1L))
  lgp <- lapply(lgp, rep, length.out = maxdim)

  border <- rep(border, length.out = maxdim)
  col <- rep(col, length.out = maxdim)
  lwd <- rep(lwd, length.out = maxdim)
  lty <- rep(lty, length.out = maxdim)

  # preallocate result
  lst <- vector("list", maxdim)

  for (i in seq_len(maxdim)) {
    # only recompute theta if nv changes
    theta <- seq(0, 2 * pi, length.out = lgp$nv[i] + 1L)[-1L]
    ct <- cos(theta)
    st <- sin(theta)

    # unrotated ellipse
    dx <- lgp$radius.x[i] * ct
    dy <- lgp$radius.y[i] * st

    # rotation
    if (lgp$rot[i] != 0) {
      cr <- cos(lgp$rot[i])
      sr <- sin(lgp$rot[i])
      ptx <- lgp$x[i] + cr * dx - sr * dy
      pty <- lgp$y[i] + sr * dx + cr * dy
    } else {
      ptx <- lgp$x[i] + dx
      pty <- lgp$y[i] + dy
    }

    # plot if requested
    if (plot) {
      polygon(ptx, pty,
        border = border[i],
        col = col[i], lty = lty[i], lwd = lwd[i]
      )
    }

    lst[[i]] <- list(x = ptx, y = pty)
  }

  # simplify result if only one ellipse
  if (maxdim == 1L) {
    lst <- grDevices::xy.coords(lst[[1]])
  } else {
    lst <- lapply(lst, grDevices::xy.coords)
  }
  invisible(lst)
}




# Modes from a kde distribution ------------------------------------------------
modes <- function(kde) {
  c(
    kde$x[kde$x < 0][which.max(kde$y[kde$x < 0])],
    kde$x[kde$x > 0][which.max(kde$y[kde$x > 0])]
  )
}


.bind_cols <- function(x, ...) {
  dots <- list(...)

  # keep only non-NULL and non-zero-length columns
  dots <- dots[vapply(dots, function(y) !is.null(y) && length(y) > 0, logical(1))]

  # convert vectors to column matrices
  dots <- lapply(dots, function(y) {
    if (is.vector(y) && !is.matrix(y)) y <- matrix(y, ncol = 1)
    y
  })

  # check row consistency
  n <- nrow(x)
  for (i in seq_along(dots)) {
    if (nrow(dots[[i]]) != n) {
      stop("All inputs must have the same number of rows as x")
    }
  }

  # combine
  as.data.frame(cbind(x, do.call(cbind, dots)), check.names = FALSE)
}



#' List of vectors
#' 
#' Creates a list of Cartesian vectors  from an spherical objects. This is a convenience 
#' function to link with the package `geologyGeometry` by J. R. Davis
#' 
#' @inheritParams sph_mean
#' 
#' @returns list
#' @export
#' 
#' @examples
#' data(example_lines)
#' l <- Line(example_lines$trend, example_lines$plunge)
#' vec_list(l) 
vec_list <- function(x){
  Vec3(x) |> 
    asplit(MARGIN = 1)
}
