# Fabric intensities and plots -------------------------------------------------

#' Fabric intensity and shape indices
#'
#' @inheritParams shape_params
#' @returns numeric vector containing the fabric shape and intensity indices:
#' \describe{
#' \item{`P`}{Point (Vollmer 1990). Range: (0, 1)}
#' \item{`G`}{Girdle (Vollmer 1990). Range: (0, 1)}
#' \item{`R`}{Random (Vollmer 1990). Range: (0, 1)}
#' \item{`B`}{Cylindricity (Vollmer 1990). Range: (0, 1)}
#' \item{`C`}{Cylindricity or Fabric strength (Woodcock 1977). Range: (0, `Inf`)}
#' \item{`I`}{Cylindricity or Fabric intensity (Lisle 1985). Range: (0, 5)}
#' \item{`D`}{"Distance" from uniformity, linear from R to P, and R to G
#' (Vollmer 2020). Range: (0, 1). End members are: uniform D = 0, girdle D = 0.5,
#' cluster D = 1. The 99% level for a test against uniformity for a sample size
#' of 300 is D = 0.1.}
#' }
#' @export
#'
#' @seealso [shape_params()]
#'
#' @examples
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' fabric_indexes(x)
fabric_indexes <- function(x) {
  shape_params(x)$Vollmer
}


#' Fabric plot of Vollmer (1990)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @inheritParams woodcock_plot
#' @param ngrid integer or 3-element vector specifying the amount of gridlines
#' for the P, G, and G axes. Constant grid spacing when only one integer is given.
#' `NULL` when no grid.
#'
#' @references Vollmer, F. W. (1990). An application of eigenvalue methods to 
#' structural domain analysis. Geological Society of America Bulletin, 102, 786<U+2013>791.
#'
#' @seealso [fabric_indexes()]
#' @family fabric-plot
#' @returns plot and when stored as an object, the `P`, `G`, and `R` values as a numeric vector.
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' a <- rvmf(100, mu = mu, k = 10)
#' vollmer_plot(a, lab = "VMF")
#'
#' set.seed(20250411)
#' b <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
#' vollmer_plot(b, lab = "FB", add = TRUE, col = "red")
#'
#' set.seed(20250411)
#' c <- runif.spherical(n = 100, "Line", method = "rotasym")
#' vollmer_plot(c, lab = "UNIF", add = TRUE, col = "green")
#'
#' set.seed(20250411)
#' d <- rkent(100, mu = mu, k = 10, b = 4)
#' vollmer_plot(d, lab = "KENT", add = TRUE, col = "blue")
#' title("Fabric plot of Vollmer (1990)")
vollmer_plot <- function(x, labels = NULL, add = FALSE, ngrid = c(5, 5, 5), ...) {
  b <- NULL
  x_vollmer <- fabric_indexes(x)
  R <- x_vollmer["R"]
  P <- x_vollmer["P"]
  G <- x_vollmer["G"]

  vec <- c(P = P, G = G, R = R)

  A <- c(0, 0) # left
  B <- c(1, 0) # right
  C <- c(1 / 2, sqrt(3) / 2) # top
  abc <- rbind(A, B, C)

  PGR <- P + G + R
  PGR <- sum(vec)

  coords <- colSums(vec * abc) / PGR

  if (!add) {
    graphics::par(xpd = TRUE)
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

  invisible(x_vollmer)
}

#' Fabric plot of Woodcock (1977)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @inheritParams ot_eigen
#' @param labels character. text labels
#' @param add logical. Should data be plotted to an existing plot?
#' @param max numeric. Maximum value for x and y axes. If `NULL`, it is calculated from the data.
#' @param ... optional graphical parameters
#'
#' @references Woodcock, N. H. (1977). Specification of fabric shapes using an eigenvalue method. Geological Society of America Bulletin88, 1231<U+2013>1236. http://pubs.geoscienceworld.org/gsa/gsabulletin/article-pdf/88/9/1231/3418366/i0016-7606-88-9-1231.pdf
#'
#' @seealso [fabric_indexes()], [ot_eigen()]
#' @family fabric-plot
#' 
#' @return A plot and when stored as an object, the orientation tensor's
#' eigenvalues and eigenvectors as a list.
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' woodcock_plot(x, lab = "x", main = "Fabric plot of Woodcock (1977)")
#' y <- rvmf(100, mu = mu, k = 20)
#' woodcock_plot(y, lab = "y", add = TRUE, col = "red")
woodcock_plot <- function(x, labels = NULL, add = FALSE, max = 7, ...) {
  x_eigen <- ot_eigen(x, scaled = TRUE)

  max_val <- if (is.null(max)) {
    max(log(x_eigen$values[1] / x_eigen$values[2]), log(x_eigen$values[2] / x_eigen$values[3])) + 1
  } else {
    max
  }

  if (!add) {
    graphics::plot(
      0, 0,
      type = "n",
      # asp = 1,
      xlab = expression(log(lambda[2] / lambda[3])),
      ylab = expression(log(lambda[1] / lambda[2])),
      xaxs = "i", yaxs = "i",
      xlim = c(0, max_val), ylim = c(0, max_val), ...
    )

    # for (j in seq(2, 8, 2)) {
    #   graphics::abline(a = j, b = -1, col = "grey", lty = 3)
    # }

    for (i in c(.2, .5, 2, 5)) {
      graphics::abline(a = 0, b = i, col = "grey80", lty = 2)
    }

    rr_breaks <- pretty(c(0, max_val))
    arc_theta <- seq(0, pi / 2, length.out = 200)
    for (rr in rr_breaks) {
      graphics::lines(rr * cos(arc_theta), rr * sin(arc_theta), col = "grey80", lty = 3)
    }

    graphics::abline(a = 0, b = 1, col = "grey", lty = 1, lwd = 1.5)
  }
  if (!is.null(labels)) {
    graphics::text(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), label = labels, ...)
  } else {
    graphics::points(log(x_eigen$values[2] / x_eigen$values[3]), log(x_eigen$values[1] / x_eigen$values[2]), ...)
  }

  invisible(x_eigen)
}

#' Fabric plot of Hsu (1965)
#'
#' @inheritParams principal_stretch
#' @inheritParams hsu_plot
#' @param ... optional parameters passed to [hsu_plot()]
#'
#' @returns plot and when stored as object, a list containing the Lode parameter `lode` and the natural octahedral strain `es`.
#' @family fabric-plot
#' @seealso [ell_lode()] for Lode parameter, and [ell_nadai] for natural octahedral strain.
#' @export
#' 
#' @references Hsu, T. C. (1966). The characteristics of coaxial and non-coaxial 
#' strain paths. Journal of Strain Analysis, 1(3), 216–222. 
#' \doi{10.1243/03093247V013216}
#'
#' @examples
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' hsu_fabric_plot(x, labels = "x")
#'
#' set.seed(20250411)
#' y <- rvmf(100, mu = mu, k = 20)
#' hsu_fabric_plot(y, labels = "y", col = "red", add = TRUE)
hsu_fabric_plot <- function(x, labels = NULL, add = FALSE, es.max = 3, ...) {
  x_eigen <- principal_stretch(x)

  # X <- x_eigen$values[1]
  # Y <- x_eigen$values[2]
  # Z <- x_eigen$values[3]
  #
  # R_YZ <- Y / Z
  # R_XY <- X / Y
  #
  # R_XZ <- R_XY * R_YZ
  # es <- 1 / sqrt(3) * sqrt(log(R_XY)^2 + log(R_YZ)^2 + log(1 / R_XZ)^2) # Nadai, 1963
  #
  #
  # K <- log(R_XY) / log(R_YZ) # Hossack 1968
  # lode <- (1 - K) / (1 + K)

  es <- ell_nadai(x_eigen)
  es.max <- if (is.null(es.max)) max(es) else es.max

  lode <- ell_lode(x_eigen)

  if (!add) {
    hsu_plot(0, 0, es.max = es.max)
  }

  # Map Lode parameter (-1..1) to angle in radians (-30°..+30° around vertical)
  theta <- lode * (pi / 6) # -1 -> -30°, 0 -> 0° (vertical), +1 -> +30°

  # Shift so plane strain = vertical (pi/2)
  theta_shift <- theta + pi / 2

  # Cartesian coordinates
  x <- es * cos(theta_shift)
  y <- es * sin(theta_shift)

  if (!is.null(labels)) {
    graphics::text(x, y, label = labels, ...)
  } else {
    graphics::points(x, y, ...)
  }

  invisible(list(lode = lode, es = es))
}
