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

rad2deg <- function(rad) {
  rad * 180 / pi
}

deg2rad <- function(deg) {
  deg * pi / 180
}

DEG2RAD <- function() {
  pi / 180
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

assign_col <- function(x, n = length(x), pal = viridis::viridis, ...) {
  normalized_data <- normalize(x)
  colors <- do.call(pal, args = list(n = n, ...))
  colors[as.numeric(cut(normalized_data, breaks = n))]
}

bin_color <- function(x, breaks, pal = viridis::viridis, ...) {
  breaks <- pretty(x, n = breaks)
  n2 <- length(breaks) - 1
  cols <- do.call(pal, args = list(n2, ...)) # [order]
  named_cols <- cut(x, breaks = breaks, labels = cols, include.lowest = TRUE) |>
    as.character()
  names(named_cols) <- cut(x, breaks = breaks, include.lowest = TRUE)
  named_cols
}



# Modes from a kde distribution ------------------------------------------------
modes <- function(kde){
  c(kde$x[kde$x < 0][which.max(kde$y[kde$x < 0])],
    kde$x[kde$x > 0][which.max(kde$y[kde$x > 0])])
}