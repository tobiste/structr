# deg2rad <- function(x) x * pi / 180
# rad2deg <- function(x) x * 180 / pi
# 
# sind <- function(x) sin(deg2rad(x))
# cosd <- function(x) cos(deg2rad(x))
# tand <- function(x) tan(deg2rad(x))
# asind <- function(x) rad2deg(asin(x))
# acosd <- function(x) rad2deg(acos(x))
# atand <- function(x) rad2deg(atan(x))
# atan2d <- function(y, x) rad2deg(atan2(y, x)) %% 360




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
#' @family parse-orientations
#' @name split
#'
#' @examples
#' test <- c("45NW", "4SE")
#' split_strike(test)
NULL

#' @rdname split
#' @export
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

#' @rdname split
#' @export
split_strike <- function(x) split_trailing_letters(x)

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
#' @param n_directions either 4, for 4-points (N, E, S, W), 8 for 8-point (N, NE, E, …) or 6 for 16-point (N, NNE, NE, …) cardinal version.
#'
#' @returns character vector
#' @export
#'
#' @family parse-orientations
#'
#' @examples
#' azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 4) # 8-point compass
#' azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 8) # 8-point compass
#' azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 16) # 16-point compass
azimuth_to_cardinal <- function(x, n_directions = 8) {
  # Normalize to 0–360
  azimuth <- x %% 360


  if (n_directions == 4) {
    dirs <- c("N", "E", "S", "W")
  } else if (n_directions == 8) {
    dirs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  } else if (n_directions == 16) {
    dirs <- c(
      "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
      "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
    )
  } else {
    stop("n_directions must be 4, 8 or 16")
  }

  step <- 360 / n_directions
  # Round to nearest sector center
  idx <- floor((azimuth + step / 2) / step) + 1
  idx[idx > n_directions] <- 1

  dirs[idx]
}


#' Converts strike and dip quadrant notation into into dip direction
#'
#' @param strike numeric. Strike in degree (left- or right-hand rule)
#' @param dip_quadrant character. Quadrant of dip direction
#' @param n_directions integer.
#'
#' @returns Dip direction in degrees
#' @family parse-orientations
#' @export
#'
#' @examples
#' s <- c(270, 315, 0, 45, 90, 135, 180, 225, 270) # strike in left-hand-rule
#' q <- c("N", "E", "E", "S", "S", "W", "W", "N", "N") # dip quadrant
#' quadrant2dd(s, q)
quadrant2dd <- function(strike, dip_quadrant, n_directions = c(4L, 8L, 16L)) {
  # 1. try right-hand-rule
  dd <- strike - 90

  # check if already correct
  res1 <- azimuth_to_cardinal(dd, n_directions = min(2 * n_directions, 16))
  match <- sapply(seq_along(res1), function(i) grepl(dip_quadrant[i], res1[i]))
  # do left-handrule if not
  ifelse(match, dd, strike + 90) %% 360
}



#' Quadrant measurement expressions to angles
#'
#' Interprets quadrant measurement expressions, such as "E30N" or "W10S" as azimuth angles
#'
#' @param x character. Dip direction and quadrant conventions
#'
#' @return numeric. Dip direction in degrees
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



# Color, shape, and size assignment helper for plotting --------------------------------------------

#' @keywords internal
#' @noRd
.normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' @keywords internal
#' @noRd
.scale <- function(x, from = range(x), to) {
  original_min <- from[1]
  original_max <- from[2]
  
  to <- sort(to)
  
  target_min <- to[1]
  target_max <- to[2]
  
  target_min + (x - original_min) * (target_max - target_min) / (original_max - original_min)
}


#' Helper functions to assign plotting colors to a vector
#'
#' @param x vector to colorize
#' @param breaks integer giving the desired number of intervals. Non-integer values are rounded down.
#' @param n integer. The number of colors (\eqn{\ge}1) to be in the palette.
#' @param title character. Legend title
#' @param pal color function; Default is [viridis::viridis()]
#' @param fill color vector
#' @param cex character expansion factor relative to current par("cex"). 
#' Used for text in legend.
#' @param legend character vector. Names of discrete colors.
#' Can be ignored when `cols` is a named vector.
#' @param position Legend position. Either a two-column vector of the x and y coordinates, or a
#' keyword from the list `"bottomright"`, `"bottom"`, `"bottomleft"`, `"left"`, 
#' `"topleft"`, `"top"`, `"topright"`, `"right"` and `"center"`.
#' @param ... arguments passed to color function
#'
#' @return character vector of colors in hexadecimal code
#' @import viridis
#' @importFrom grDevices colorRamp
#' @seealso [PlotTools::SpectrumLegend()] - salternative tool to generate a 
#' color-bar in base-R plots
#' @family assign
#' @name assign-color
#'
#' @examples
#' set.seed(20250411)
#'
#' # example for discrete colors
#' x <- rvmf(5, mu = Line(120, 50), k = 5)
#' key <- letters[round(runif(5, 1, 26))]
#' plot(x, col = assign_col_d(key), grid.params = list(guides = FALSE))
#' legend_col_d(assign_col_d(key))
#'
#' # example for continuous colors:
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' plot(x, col = assign_col(runif(100)), grid.params = list(guides = FALSE))
#' legend_col(seq(0, 1, .1), title = "test")
NULL

#' @rdname assign-color
#' @export
assign_col_d <- function(x, pal = viridis::viridis, ...) {
  groups <- unique(x)
  n <- length(groups)
  cols <- do.call(pal, args = list(n = n, ...))
  named_cols <- stats::setNames(cols, groups)
  named_cols[x]
}

#' @rdname assign-color
#' @export
assign_col <- function(x, n = length(x), pal = viridis::viridis, ...) {
  normalized_data <- .normalize(x)
  colors <- do.call(pal, args = list(n = n, ...))
  colors[as.numeric(cut(normalized_data, breaks = n))]
}

#' @rdname assign-color
#' @export
assign_col_binned <- function(x, breaks = 5, pal = viridis::viridis, ...) {
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

#' @rdname assign-color
#' @export
legend_col <- function(breaks, title = NULL, pal = viridis::viridis, cex = 1, ...) {
  label_pos <- .normalize(breaks)
  legend_image <- grDevices::as.raster(
    matrix(
      color_func(seq(0, 1, 0.1), pal = pal, ...),
      ncol = 1
    )
  )

  graphics::par(new = TRUE, xpd = NA)
  graphics::layout(matrix(1, 1))

  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
  graphics::text(x = 0.9 + (0.1 * 1.1), y = label_pos, labels = round(breaks, 2), adj = 0, cex = .8 * cex * par("cex"))
  graphics::text(x = .925, y = .5, labels = title, adj = 0.5, srt = 90, font = 2, cex = cex * par("cex"))
  graphics::rasterImage(legend_image, .95, 1, 1, 0)
}

#' @rdname assign-color
#' @export
legend_col_d <- function(fill, legend = names(fill), position = "topright", ...) {
  if(length(position)==1){
    xpos <- position
    ypos <- NULL
  } else {
    xpos <- position[1]
    ypos <- position[2]
  }
  
  graphics::legend(position,
    legend = legend,
    fill = fill,
    ...
  )
}



#' Assign plotting size (cex values) to a vector
#' 
#' `assign_cex()` maps the character expansion (size). 
#'  The cex is most commonly used for points and text, and humans perceive the 
#'  area of points (not their radius), so this provides for optimal perception. 
#'  The argument `area` ensures that a value of `0` is mapped to a size of `0`; 
#' `assign_cex_binned()` is a binned version, and `assign_cex_d()` assigns 
#' cex values to discrete values.
#' 
#' @details The character expansion `cex` is a number indicating the amount by 
#' which plotting text and symbols 
#' should be scaled relative to the default. `1`=default, `1.5` is 50% larger, 
#' `0.5` is 50% smaller, etc.
#'
#' @param x vector
#' @param range numeric 2-element vector. Output range of `cex` values. Minimum 
#' value must be greater than 0.
#' @param pch plotting character to be used in legend.
#' @param area logical. Whether `cex` should be proportional to the area (`TRUE`) 
#' or the radius (`FALSE`, the default) of the plotting character.
#' @param values numeric. `cex` values to manually assign to `x`. Must be at least 
#' the number of unique values in `x`.
#' @inheritParams assign_col
#' @param ... arguments passed to [graphics::legend()]
#' 
#' @family assign
#' @seealso [PlotTools::SizeLegend()] - alternative tool to generate nice 
#' looking legends for cex values
#'
#' @returns numeric vector
#' @name assign-cex
#'
#' @examples
#' set.seed(20250411)
#'
#' # example for continuous colors:
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' key <- runif(100)
#' plot(x, cex = assign_cex(key), grid.params = list(guides = FALSE))
#' legend_cex(key, position = 'topright', area = TRUE)
NULL

#' @rdname assign-cex
#' @export
assign_cex <- function(x, range = c(0.25, 2), area = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(min(range) > 0)
  
  if(isTRUE(area)) x <- .scale(sqrt(abs(x)), c(0, max(x)), c(0, 1))
  .scale(x, to = range)
}

#' @rdname assign-cex
#' @export
assign_cex_binned <- function(x, range = c(0.25, 2), breaks = 5, area = FALSE){
  x_breaks <- pretty(x, n = breaks)
  n2 <- length(breaks) - 1
  cexs <- assign_cex(seq_len(n2), range, area)
  named_cexs <- cut(x, breaks = x_breaks, labels = cexs, include.lowest = TRUE) 
  names(named_cexs) <- cut(x, breaks = x_breaks, include.lowest = TRUE)
  named_cexs
}

#' @rdname assign-cex
#' @export
assign_cex_d <- function(x, values = NULL, range = c(0.25, 2)){
  xf <- as.factor(x)
  n <- nlevels(xf)
  
  range <- sort(range)
  
  # If user provides custom values:
  if (!is.null(values)) {
    if (length(values) < n) {
      stop("Length of 'values' must be at least the number of unique values in 'x'.")
    }
    # Use first n values
    cex_vals <- values[seq_len(n)]
  } else {
    # Otherwise generate cex values across default_range
    cex_vals <- seq(range[1], range[2], length.out = n)
  }
  
  # Assign to each element of x based on factor indexing
  assigned <- cex_vals[as.integer(xf)]
  names(assigned) <- values
  return(assigned)
}

#' @rdname assign-cex
#' @export
legend_cex <- function(x, range = c(0.25, 2), breaks = 5, values = NULL, area = FALSE, position = "topright", pch = 16, ...){
    if(!is.null(values)) {
      cexs <- assign_cex_d(x, values, range)
    } else {
      x_breaks <- pretty(x, n = breaks)
      n2 <- length(x_breaks) - 1
      cexs <- assign_cex(seq_len(n2), range, area)
      names(cexs) <- cut(x, breaks = x_breaks, include.lowest = TRUE, ordered_result = TRUE) |> 
       levels()
    }
  
  if(length(position)==1){
    xpos <- position
    ypos <- NULL
  } else {
    xpos <- position[1]
    ypos <- position[2]
    }
  
  graphics::legend(x = xpos, y = xpos,
                   legend = unique(names(cexs)),
                   pt.cex = unique(cexs),
                   pch = pch,
                   ...
  )
}


#' Assigns plotting characters (pch values) to a vector
#' 
#' `assign_pch()` maps discrete variables to six easily discernible shapes. If you have more 
#' than six levels, you will get a warning message, and the seventh and 
#' subsequent levels will not appear on the plot. 
#' You can not map a continuous variable to shape unless `assign_pch_binned()` is 
#' used.
#'
#' @param x vector.
#' @param solid Should the plotting character be solid, `TRUE` (the default), or hollow, `FALSE`?
#' @inheritParams assign_col
#' @param ... arguments passed to `graphics::legend()`
#' 
#' @family assign
#'
#' @returns named integer vector
#' @name assign-pch
#'
#' @examples
#' set.seed(20250411)
#'
#' # example for discrete colors
#' x <- rvmf(5, mu = Line(120, 50), k = 5)
#' key <- sample(letters, 5, replace = TRUE)
#' plot(x, pch = assign_pch(key), grid.params = list(guides = FALSE))
#' legend_pch(key)
NULL

#' @rdname assign-pch
#' @export
assign_pch <- function(x, solid = TRUE){
  x_unique <- unique(x) |> sort()
  xn <- length(x_unique) 
  if(isTRUE(solid)){
    shapes <- c(16, 17, 15, 3, 7, 8)
  } else {
    shapes <- c(1, 2, 0, 3, 7, 8)
  }
  
  if(length(x_unique) > length(shapes)) warning(
    paste0(
      "The shape palette can deal with a maximum of 6 discrete values because more than 6 becomes difficult to discriminate.\n You have requested ", xn, " values. Consider specifying shapes manually if you need that many of them.")
    )
  
  xf <- as.factor(x)
  assigned <- shapes[(as.integer(xf) - 1) %% length(shapes) + 1]
  names(assigned) <- xf
  
  return(assigned)
}

#' @rdname assign-pch
#' @export
assign_pch_binned <- function(x, solid = TRUE, breaks = 6){
  breaks <- pretty(x, n = breaks)
  n2 <- length(breaks) - 1
  
  pchs <- assign_pch(seq_len(n2), solid)
  
  named_pchs <- cut(x, breaks = breaks, labels = pchs, include.lowest = TRUE) |>
    as.character()
  names(named_pchs) <- cut(x, breaks = breaks, include.lowest = TRUE)
  named_pchs
}

#' @rdname assign-pch
#' @export
legend_pch <- function(x, solid = TRUE, breaks = NULL, position = "topright", ...){
  pchs <- if(is.null(breaks)) assign_pch(x, solid) else assign_pch_binned(x, solid, breaks)
  
  if(length(position)==1){
    xpos <- position
    ypos <- NULL
  } else {
    xpos <- position[1]
    ypos <- position[2]
  }
  
  graphics::legend(position,
                   legend = unique(names(pchs)),
                   pch = unique(pchs),
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
#' @keywords internal
#' @noRd
modes <- function(kde) {
  c(
    kde$x[kde$x < 0][which.max(kde$y[kde$x < 0])],
    kde$x[kde$x > 0][which.max(kde$y[kde$x > 0])]
  )
}

#' @keywords internal
#' @noRd
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
#' @param ls list of 3-element vectors (Cartesian coordinates x, y, z)
#'
#' @returns `vec_list` returns a list of 3-element vectors (Cartesian coordinates x, y, z).
#' `list_vec` returns a `"Vec3"` object
#'
#' @name vec-list
# #' @seealso [rot2pair()] and [pair2rot()]
#'
#' @examples
#' ls <- vec_list(example_lines[1:5, ])
#' print(ls)
#'
#' list_vec(ls)
NULL

#' @rdname vec-list
#' @export
vec_list <- function(x) {
  Vec3(x) |>
    asplit(MARGIN = 1)
}

#' @rdname vec-list
#' @export
list_vec <- function(ls) {
  do.call(rbind, ls) |>
    as.Vec3()
}

# Confidence margin for scalar values
# .confidence_margin_scalar <- function(x, alpha = 0.05){
#   p <- 1 - alpha/2
#   stats::qt(p,df=n-1)*s/sqrt(n)
# }


#' List transposition
#'
#' Rearranges a list of n matrices with n rows and 3 columns into a list of n
#' matrices with n rows and three columns.
#'
#' @param x list of 3-column matrices with each same number of rows
#'
#' @return list
#' @noRd
#'
#' @examples
#' d <- defgrad_from_generalshear(k = 2.5, gamma = 0.9)
#' v <- velgrad(d, time = 10)
#' d_steps <- defgrad(v, time = 10, steps = 2)
#'
#' # apply on orientation data
#' set.seed(20250411)
#' l <- rvmf(100, mu = Line(0, 90), k = 100)
#' l_trans <- lapply(d_steps, function(i) {
#'   transform_linear(l, i)
#' })
#'
#' transpose_list(l_trans)
transpose_list <- function(x) {
  # Check that input is a non-empty list
  stopifnot(is.list(x), length(x) > 0)

  # Extract all dimensions
  dims <- lapply(x, dim)

  # Check all dimensions are identical
  stopifnot(all(vapply(dims, function(d) all(d == dims[[1]]), logical(1))))

  # Extract i (rows) and n (number of matrices)
  i <- dims[[1]][1]
  n <- length(x)

  # Transpose list: create i matrices, each with n rows and 3 columns
  x_rearranged <- lapply(seq_len(i), function(r) {
    do.call(rbind, lapply(x, function(M) {
      M <- unclass(M)
      M[r, , drop = FALSE]
    }))
  })

  return(x_rearranged)
}
