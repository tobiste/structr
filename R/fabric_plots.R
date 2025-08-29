# Fabric intensities and plots -------------------------------------------------

#' Fabric intensity and shape indices
#'
#' @inheritParams best_fit_plane
#' @returns numeric vector containing the fabric shape and intensity indices:
#' \describe{
#' \item{`P`}{Point (Vollmer 1990). Range: (0, 1)}
#' \item{`G`}{Girdle (Vollmer 1990). Range: (0, 1)}
#' \item{`R`}{Random (Vollmer 1990). Range: (0, 1)}
#' \item{`B`}{cylindricity (Vollmer 1990). Range: (0, 1)}
#' \item{`C`}{cylindricity or Fabric strength (Woodcock 1977). Rrange: (0, Inf)}
#' \item{`I`}{cylindricity or Fabric intensity (Lisle 1985). Range: (0, 5)}
#' \item{`D`}{"distance" from uniformity, linear from R to P, and R to G (Vollmer 2020). Range: (0, 1). End members are: uniform D = 0, girdle D = 0.5, cluster D = 1. The 99% level for a test against uniformity for a sample size of 300 is D = 0.1.}
#' }
#' @export
#'
#' @seealso [or_shape_params()]
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' fabric_indexes(x)
fabric_indexes <- function(x) {
  or_shape_params(x)$Vollmer
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
#' c <- runif.spherical(n = 100, "Line", method = "rotasym")
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

#' Fabric plot of Woodcock (1977)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @inheritParams best_fit_plane
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
  x_eigen <- eigen(x, scaled = TRUE)
  if (!is.null(labels)) {
    graphics::text(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), label = labels, ...)
  } else {
    graphics::points(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), ...)
  }
}
