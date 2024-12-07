# Create a ternary plot from a data frame
ternaryPlot <- function(x, plotPoints = TRUE, labels = c("", "", ""), grid = TRUE, increment = 20, ...) {
  ternaryTriangle()
  ternaryLabels(labels[1], labels[2], labels[3])
  if (grid == TRUE) {
    ternaryGrid(increment)
  }
  if (plotPoints == TRUE) {
    ternaryPoints(x, ...)
  }
}

# Add points from a data frame to an existing ternary plot
ternaryPoints <- function(x, ...) {
  x <- validatedTernaryPoints(x)
  coords <- cartesianFromTernary(x[, 1], x[, 2], x[, 3])
  graphics::points(coords$x, coords$y, ...)
}

# Add a line segment to an existing ternary plot
# Color and line type may be added as additional arguments
ternarySegment <- function(x0, x1, ...) {
  # x0 and x1 are vectors of the endpoint ternary coordinates
  coords0 <- cartesianFromTernary(x0[1], x0[2], x0[3])
  coords1 <- cartesianFromTernary(x1[1], x1[2], x1[3])
  graphics::segments(coords0$x, coords0$y, coords1$x, coords1$y, ...)
}

# Add a polygon to an existing ternary plot
# Color may be added as an additional argument
ternaryPolygon <- function(x, ...) {
  nPoints <- nrow(x)
  xCoord <- vector(mode = "numeric", length = nPoints)
  yCoord <- vector(mode = "numeric", length = nPoints)
  for (i in 1:nPoints) {
    coords <- cartesianFromTernary(x[i, 1], x[i, 2], x[i, 3])
    xCoord[i] <- coords$x
    yCoord[i] <- coords$y
  }
  graphics::polygon(xCoord, yCoord, ...)
}

# Add text to an existing ternary plot
# Text styling may be added as additional arguments
ternaryText <- function(x, label = "", ...) {
  coords <- cartesianFromTernary(x[1], x[2], x[3])
  graphics::text(coords$x, coords$y, label = label, ...)
}

# ---------------------------------------------------------------------------------------
# The following functions are called by ternaryPlot() and generally will not need to be
# called directly

## Plotting primitives -------------------------------------------------------------------
ternaryTriangle <- function() {
  top <- cartesianFromTernary(100, 0, 0)
  left <- cartesianFromTernary(0, 100, 0)
  right <- cartesianFromTernary(0, 0, 100)
  lim <- c(-1.1, 1.1)
  graphics::plot(top$x, top$y, xlim = lim, ylim = lim, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "")
  graphics::segments(top$x, top$y, right$x, right$y)
  graphics::segments(top$x, top$y, left$x, left$y)
  graphics::segments(left$x, left$y, right$x, right$y)
}

ternaryLabels <- function(top = "", left = "", right = "") {
  topCoord <- cartesianFromTernary(100, 0, 0)
  leftCoord <- cartesianFromTernary(0, 100, 0)
  rightCoord <- cartesianFromTernary(0, 0, 100)
  graphics::text(topCoord$x, topCoord$y, top, pos = 3)
  graphics::text(leftCoord$x, leftCoord$y, left, pos = 1, srt = -60)
  graphics::text(rightCoord$x, rightCoord$y, right, pos = 1, srt = 60)
}

ternaryGrid <- function(increment) {
  low <- increment
  high <- 100 - increment

  m <- seq(low, high, increment)
  nLines <- length(m)

  n1 <- o2 <- seq(high, low, -increment)
  n2 <- o1 <- rep(0, nLines)

  for (i in 1:nLines) {
    a <- cartesianFromTernary(m[i], n1[i], o1[i])
    b <- cartesianFromTernary(m[i], n2[i], o2[i])
    graphics::segments(a$x, a$y, b$x, b$y, col = "lightgray", lty = 3)

    a <- cartesianFromTernary(n1[i], m[i], o1[i])
    b <- cartesianFromTernary(n2[i], m[i], o2[i])
    graphics::segments(a$x, a$y, b$x, b$y, col = "lightgray", lty = 3)

    a <- cartesianFromTernary(n1[i], o1[i], m[i])
    b <- cartesianFromTernary(n2[i], o2[i], m[i])
    graphics::segments(a$x, a$y, b$x, b$y, col = "lightgray", lty = 3)
  }
}

## Convert from ternary coordinates to cartesian (x, y) coordinates ----------------------
cartesianFromTernary <- function(top, left, right) {
  y <- (top - 50) / 50 # vertically spans from -1 to 1
  baseHalfWidth <- 2 / tan(60 * pi / 180) # : equilateral triangle
  horizontalHalfWidth <- ((100 - top) * baseHalfWidth) / 100
  horizontalProportion <- (right / (right + left + 0.0000001) - 0.5) * 2
  x <- horizontalProportion * horizontalHalfWidth
  xyCoords <- data.frame(cbind(x = x, y = y))
  colnames(xyCoords) <- c("x", "y")
  xyCoords
}

## Remove rows with NA values and rows that don't sum to 100 -----------------------------
validatedTernaryPoints <- function(x) {
  rowsWithNA <- rowSums(is.na(x))
  xNAremoved <- x[rowsWithNA == 0, ]
  rowsEqualing100 <- abs(rowSums(xNAremoved) - 100) <= 2
  xFinal <- xNAremoved[rowsEqualing100, ]
  xFinal
}

# Fabric intensities and plots -------------------------------------------------


#' Fabric intensity and shape indices
#'
#' @param x numeric. Can be three element vector, three column array, or an object of class `"line"` or `"plane"`
#'
#' @returns numeric vector containing the fabric shape and intensity indices:
#' \describe{
#' \item{`P`}{Point (Vollmer 1990), range: (0, 1)}
#' \item{`G`}{Girdle (Vollmer 1990), range: (0, 1)}
#' \item{`R`}{Random (Vollmer 1990), range: (0, 1)}
#' \item{`B`}{cylindricity (Vollmer 1990), range: (0, 1)}
#' \item{`C`}{cylindricity or Fabric strength (Woodcock 1977), range: (0, Inf)}
#' \item{`I`}{cylindricity or Fabric intensity (Lisle 1985), range: (0, 5)}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' fabric_indexes(x)
fabric_indexes <- function(x) {
  x_eigen <- or_eigen(x)

  N <- nrow(x)
  P <- (x_eigen$values[1] - x_eigen$values[2]) / N * 100
  G <- 2 * (x_eigen$values[2] - x_eigen$values[3]) / N * 100
  R <- 3 * (x_eigen$values[3]) / N * 100

  B <- (P + G)
  C <- log(x_eigen$values[1] / x_eigen$values[3])
  I <- 7.5 * ((x_eigen$values[1] / N - (1 / 3))^2 + (x_eigen$values[2] / N - (1 / 3))^2 + (x_eigen$values[3] / N - (1 / 3))^2)

  c(P = P, G = G, R = R, B = B, C = C, I = I)
}


#' Fabric plot of Vollmer (1990)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @inheritParams WoodcockPlot
#' @param ngrid integer or 3-element vector specifying the amount of gridlines
#' for the P, G, and G axes. Constant grid spacing when only one integer is given.
#' `NULL` when no grid.
#'
#' @references Vollmer, F. W. (1990). An application of eigenvalue methods to structural domain analysis. Geological Society of America Bulletin, 102, 786<U+2013>791.
#'
#' @seealso [WoodcockPlot()], [fabric_indexes()]
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' a <- rvmf(100, mu = mu, k = 10)
#' VollmerPlot(a, lab = "VMF")
#'
#' set.seed(1)
#' b <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
#' VollmerPlot(b, lab = "FB", add = TRUE, col = "red")
#'
#' set.seed(1)
#' c <- v_unif("line", n = 100, method = "rotasym")
#' VollmerPlot(c, lab = "UNIF", add = TRUE, col = "green")
#'
#' set.seed(1)
#' d <- rkent(100, mu = mu, k = 10, b = 4)
#' VollmerPlot(d, lab = "KENT", add = TRUE, col = "blue")
#' title("Fabric plot of Vollmer (1990)")
VollmerPlot <- function(x, labels = NULL, add = FALSE, ngrid = c(5, 5, 5), ...) {
  x_vollmer <- fabric_indexes(x)
  R <- x_vollmer["R"]
  P <- x_vollmer["P"]
  G <- x_vollmer["G"]

  vec <- c(P, G, R)

  A <- c(0, 0) # left
  B <- c(1, 0) # right
  C <- c(1 / 2, sqrt(3) / 2) # top
  abc <- rbind(A, B, C)

  PGR <- P + G + R
  PGR <- sum(vec)

  coords <- colSums(vec * abc) / PGR

  if (!add) {
    graphics::plot(1, "n", ylim = c(0, sqrt(3) / 2), xlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "")

    if (!is.null(ngrid)) {
      ngrid <- round(ngrid) + 1
      if (length(ngrid) == 1) ngrid <- rep(ngrid, 3)
      # horizontal lines
      bottom <- seq(0, 1, length.out = ngrid[1])
      bottom <- bottom[2:(length(bottom) - 1)]

      l1 <- rbind(bottom, rep(0, length(bottom))) |> t()
      rot <- cbind(c(1 / 2, sind(60)), c(-sind(60), 1 / 2))
      l2 <- t(rot %*% t(l1))
      l3 <- l2
      l3[, 1] <- 1 - l2[, 1]

      graphics::segments(l2[, 1], l2[, 2], l3[, 1], l3[, 2], col = "lightgray", lty = 3)


      # diagonal lines of left
      bottom2 <- seq(0, 1, length.out = ngrid[2])
      bottom2 <- bottom2[2:(length(bottom2) - 1)]
      l1_2 <- rbind(bottom2, rep(0, length(bottom2))) |> t()
      rot <- cbind(c(1 / 2, sind(60)), c(-sind(60), 1 / 2))
      l2_2 <- t(rot %*% t(l1_2))
      graphics::segments(l1_2[, 1], l1_2[, 2], l2_2[, 1], l2_2[, 2], col = "lightgray", lty = 3)

      # diagonal lines of right
      l2_3 <- l3
      l2_3[, 1] <- rev(l3[, 1])
      l2_3[, 2] <- rev(l3[, 2])

      graphics::segments(l1_2[, 1], l1_2[, 2], l2_3[, 1], l2_3[, 2], col = "lightgray", lty = 3)
    }

    graphics::polygon(abc)

    l1 <- colSums(c(0, 1, 1) * abc) / 2
    l2 <- colSums(c(0, 0, 1) * abc) / 2
    l3 <- c(.5, 0)

    # add axes labels
    graphics::text(l3[1], l3[2], labels = expression(lambda[3] == 0), pos = 3, offset = -1, col = "grey")
    graphics::text(l2[1], l2[2], labels = expression(lambda[2] == lambda[3]), pos = 3, srt = 60, offset = 1, col = "grey")
    graphics::text(l1[1], l1[2], labels = expression(lambda[1] == lambda[2]), pos = 3, srt = -60, offset = 1, col = "grey")

    abc_text <- abc #+ c(-.02, .02, .02)
    for (t in 1:3) graphics::text(abc_text[t, 1], abc_text[t, 2], labels = c("Point", "Girdle", "Random")[t], pos = c(1, 1, 3)[t], srt = c(-60, 60, 0)[t])
  }

  if (!is.null(labels)) {
    graphics::text(coords[1], coords[2], labels = labels, ...)
  } else {
    graphics::points(coords[1], coords[2], ...)
  }
}

# VollmerPlot <- function(x, lab = NULL, add = FALSE, grid = TRUE, increment = 20, ...) {
#   if (!add) {
#     ternaryPlot(plotPoints = FALSE, labels = c("Random", "Point", "Girdle"), grid = grid, increment = increment)
#     bottomCoord <- cartesianFromTernary(0, 50, 50)
#     leftCoord <- cartesianFromTernary(50, 50, 0)
#     rightCoord <- cartesianFromTernary(50, 0, 50)
#     graphics::text(bottomCoord$x, bottomCoord$y, expression(lambda[3] == 0), pos = 3, offset = -1, col = "grey")
#     graphics::text(leftCoord$x, leftCoord$y, expression(lambda[2] == lambda[3]), pos = 3, srt = 60, offset = 1, col = "grey")
#     graphics::text(rightCoord$x, rightCoord$y, expression(lambda[1] == lambda[2]), pos = 3, srt = -60, offset = 1, col = "grey")
#   }
#
#   x_vollmer <- fabric_indexes(x)
#   dat <- c(
#     top = x_vollmer["R"], left = x_vollmer["P"], right = x_vollmer["G"]
#   ) * 100
#
#   if (!is.null(lab)) {
#     ternaryText(dat, label = lab, cex = .8, ...)
#   } else {
#     coords <- cartesianFromTernary(dat[1], dat[2], dat[3])
#     graphics::points(coords$x, coords$y, ...)
#   }
# }


#' Fabric plot of Woodcock (1977)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @param x numeric. Can be three element vector, three column array, or an object of class `"line"` or `"plane"`
#' @param labels character. text labels
#' @param add logical. Should data be plotted to an existing plot?
#' @param ... optional graphical parameters
#' @references Woodcock, N. H. (1977). Specification of fabric shapes using an eigenvalue method. Geological Society of America Bulletin88, 1231<U+2013>1236. http://pubs.geoscienceworld.org/gsa/gsabulletin/article-pdf/88/9/1231/3418366/i0016-7606-88-9-1231.pdf
#'
#' @seealso [VollmerPlot()], [fabric_indexes()]
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' WoodcockPlot(x, lab = "x", main = "Fabric plot of Woodcock (1977)")
#' y <- rvmf(100, mu = mu, k = 20)
#' WoodcockPlot(y, lab = "y", add = TRUE, col = "red")
WoodcockPlot <- function(x, labels = NULL, add = FALSE, ...) {
  if (!add) {
    graphics::plot(
      x = 1, type = "n",
      xlab = expression(log(lambda[2] / lambda[3])),
      ylab = expression(log(lambda[1] / lambda[2])),
      xaxs = "i", yaxs = "i",
      xlim = c(0, 7), ylim = c(0, 7), ...
    )

    for (j in seq(2, 8, 2)) {
      graphics::abline(a = j, b = -1, col = "grey", lty = 3)
    }

    for (i in c(.2, .5, 2, 5)) {
      graphics::abline(a = 0, b = i, col = "grey", lty = 2)
    }
    graphics::abline(a = 0, b = 1, col = "grey", lty = 1, lwd = 1.5)
  }
  x_eigen <- or_eigen(x, scaled = TRUE)
  if (!is.null(labels)) {
    graphics::text(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), label = labels, ...)
  } else {
    graphics::points(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), ...)
  }
}
