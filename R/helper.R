vec2mat <- function(x) {
  class <- class(x)
  if (is.null(dim(x))) {
    m <- as.matrix(t(x))
  } else {
    m <- as.matrix(x)
  }
  class(m) <- class
  m
}

DEG2RAD <- function() {
  pi / 180
}

rad2deg <- function(rad) {
  rad / DEG2RAD()
}

deg2rad <- function(deg) {
  deg * DEG2RAD()
}



sind <- function(x) {
  sinpi(x / 180)
}

cosd <- function(x) {
  cospi(x / 180)
}

tand <- function(x) {
  sinpi(x / 180) / cospi(x / 180)
}


asind <- function(x) {
  asin(x) * 180 / pi
}

acosd <- function(x) {
  acos(x) * 180 / pi
}

atand <- function(x) {
  atan(x) * 180 / pi
}

atan2d <- function(x1, x2) {
  atan2(x1, x2) * 180 / pi
}

cot <- function(x) {
  1 / tan(x)
}

cotd <- function(x) {
  1 / tand(x)
}


# Color assignment helper functions --------------------------------------------

normalize <- function(x) {
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
#' stereoplot(guides = FALSE)
#' stereo_point(x, col = assign_col_d(key))
#' legend_d(assign_col_d(key))
#'
#' # example for continuous colors:
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' stereoplot(guides = FALSE)
#' stereo_point(x, col = assign_col(runif(100)))
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
  normalized_data <- normalize(x)
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
  label_pos <- normalize(breaks)
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

# Modes from a kde distribution ------------------------------------------------
modes <- function(kde) {
  c(
    kde$x[kde$x < 0][which.max(kde$y[kde$x < 0])],
    kde$x[kde$x > 0][which.max(kde$y[kde$x > 0])]
  )
}
