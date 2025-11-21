# Fabric intensities and plots -------------------------------------------------

#' Fabric Intensity and Shape of Orientation Tensor
#'
#' Fabric intensity and shape parameters of the orientation tensor based on Vollmer (1990)
#'
#' @inheritParams geodesic_mean
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
#' \item{`U`}{Uniformity statistic of Mardia (1972)}
#' }
#'
#' @references
#' Lisle, Richard J.  (1985): "The use of the orientation tensor for the description and statistical testing of fabrics." Journal of Structural Geology 7.1: 115-117.
#'
#' Mardia, Kantilal Varichand (1975): "Statistics of directional data." Journal of the Royal Statistical Society Series B: Statistical Methodology 37.3: 349-371.
#'
#' Vollmer, Frederick W. (1990): "An application of eigenvalue methods to structural domain analysis." Geological Society of America Bulletin 102.6: 786-791.
#'
#' Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on a Triangular Plot Using a Fabric Density Index and Crystal Axes Eigenvector Barycenters." Geological Society of America Abstracts. Vol. 52.
#'
#' Woodcock, Nigel H.  (1977): "Specification of fabric shapes using an eigenvalue method." Geological Society of America Bulletin 88.9: 1231-1236.
#'
#' @export
#'
#' @seealso [shape_params()], [ortensor()], [vollmer_plot()]
#'
#' @examples
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' vollmer(x)
#'
#' # Pair objects:
#' vollmer(simongomez)
vollmer <- function(x) {
  stopifnot(is.spherical(x))
  ot <- ortensor(x)
  eig <- eigen(ot, only.values = TRUE)$values #|> sort(decreasing = TRUE)

  if (is.Pair(x) & any(eig < 0)) eig <- abs(eig) |> sort(decreasing = TRUE)


  # Vollmer
  N <- nrow(x)
  P <- eig[1] - eig[2] #  Point index (Vollmer, 1990)
  G <- 2 * (eig[2] - eig[3]) #  Girdle index (Vollmer, 1990)
  R <- 3 * eig[3] # Random index (Vollmer, 1990)
  B <- P + G #  Cylindricity index (Vollmer, 1990)
  C <- log(eig[1] / eig[3])
  I <- 7.5 * sum((eig / N - 1 / 3)^2)

  U <- (15 * N / 2) * sum((eig - 1 / 3)^2) # Mardia uniformity statistic
  D <- sqrt(U / (5 * N)) # D of Vollmer 2020

  c(P = P, G = G, R = R, B = B, C = C, I = I, D = D, U = U)
}

#' Fabric Plot of Vollmer (1990)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @param x spherical object or a three-column matrix, where the first column is
#' P, the second is G, and the third one is R of the Vollmer parameters.
#' @inheritParams woodcock_plot
#' @param ngrid integer or 3-element vector specifying the amount of gridlines
#' for the P, G, and G axes. Constant grid spacing when only one integer is given.
#' `NULL` when no grid.
#'
#' @references Vollmer, F. W. (1990). An application of eigenvalue methods to
#' structural domain analysis. Geological Society of America Bulletin, 102, 786<U+2013>791.
#'
#' @seealso [vollmer()], [ortensor()]
#' @family fabric-plot
#' @returns plot and when stored as an object, the `P`, `G`, and `R` values as a numeric vector.
#' @name vollmer-plot
#'
#' @examples
#' # Orientation data
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' a <- rvmf(10, mu = mu, k = 10)
#' vollmer_plot(a, labels = "VMF")
#'
#' set.seed(20250411)
#' b <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
#' vollmer_plot(b, labels = "FB", add = TRUE, col = "red")
#'
#' set.seed(20250411)
#' c <- runif.spherical(n = 100, "Line", method = "rotasym")
#' vollmer_plot(c, labels = "UNIF", add = TRUE, col = "green")
#'
#' set.seed(20250411)
#' d <- rkent(100, mu = mu, k = 10, b = 4)
#' vollmer_plot(d, labels = "KENT", add = TRUE, col = "blue")
NULL

#' @rdname vollmer-plot
#' @export
vollmer_plot <- function(x, labels, add, ngrid, main, ...) UseMethod("vollmer_plot")

#' @rdname vollmer-plot
#' @export
vollmer_plot.default <- function(x, labels = NULL, add = FALSE, ngrid = c(5, 5, 5), main = "Vollmer diagram", ...) {
  A <- c(0, 0) # left
  B <- c(1, 0) # right
  C <- c(1 / 2, sqrt(3) / 2) # top
  abc <- rbind(A, B, C)
  coords <- sapply(seq_len(nrow(x)), function(i) {
    vec <- c(P = x[i, 1], G = x[i, 2], R = x[i, 3])
    PGR <- x[i, 1] + x[i, 3] + x[i, 3]
    PGR <- sum(vec)
    colSums(vec * abc) / PGR
  }) |> t()

  if (isFALSE(add)) .vollmer_plot_blank(ngrid, main = main)

  if (!is.null(labels)) {
    graphics::text(coords[, 1], coords[, 2], labels = labels, ...)
  } else {
    graphics::points(coords[, 1], coords[, 2], ...)
  }
}

#' @rdname vollmer-plot
#' @export
vollmer_plot.spherical <- function(x, labels = NULL, add = FALSE, ngrid = c(5, 5, 5), main = "Vollmer diagram", ...) {
  x_vollmer <- vollmer(x)
  # P <- x_vollmer["P"]
  # G <- x_vollmer["G"]
  # R <- x_vollmer["R"]
  vollmer_plot.default(
    t(x_vollmer[c("P", "G", "R")]),
    labels = labels, add = add, ngrid = ngrid, main = main, ...
  )
  invisible(x_vollmer)
}

#' @rdname vollmer-plot
#' @export
vollmer_plot.list <- function(x, labels = NULL, add = FALSE, ngrid = c(5, 5, 5), main = "Vollmer diagram", ...) {
  x_vollmer <- do.call(rbind, lapply(x, vollmer))

  vollmer_plot.default(
    x_vollmer[, c("P", "G", "R")],
    labels = labels, add = add, ngrid = ngrid, main = main, ...
  )
  invisible(x_vollmer)
}

.vollmer_plot_blank <- function(ngrid = c(5, 5, 5), main = "Vollmer diagram") {
  A <- c(0, 0) # left
  B <- c(1, 0) # right
  C <- c(1 / 2, sqrt(3) / 2) # top
  abc <- rbind(A, B, C)

  graphics::par(xpd = TRUE)
  graphics::plot(c(0, 1), c(0, sqrt(3) / 2), "n",
    asp = 1, axes = FALSE,
    main = main,
    xlab = "", ylab = ""
  )

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


#' Fabric Plot of Woodcock (1977)
#'
#' Creates a fabric plot using the eigenvalue method
#'
#' @inheritParams ot_eigen
#' @param main character. The main title for the plot.
#' @param labels character. text labels
#' @param add logical. Should data be plotted to an existing plot?
#' @param max numeric. Maximum value for x and y axes. If `NULL`, it is calculated from the data.
#' @param ... optional graphical parameters
#'
#' @references Woodcock, N. H. (1977). Specification of fabric shapes using an eigenvalue method. Geological Society of America Bulletin88, 1231<U+2013>1236. http://pubs.geoscienceworld.org/gsa/gsabulletin/article-pdf/88/9/1231/3418366/i0016-7606-88-9-1231.pdf
#'
#' @seealso [vollmer()], [ot_eigen()]
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
#' woodcock_plot(x, lab = "x")
#' y <- rvmf(100, mu = mu, k = 20)
#' woodcock_plot(y, lab = "y", add = TRUE, col = "red")
woodcock_plot <- function(x, labels = NULL, add = FALSE, max = 7, main = "Woodcock diagram", ...) {
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
      main = main,
      xpd = FALSE,
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

#' Fabric Plot of Hsu (1965)
#'
#' 3D strain diagram using the Hsu (1965) method to display the amount of the natural
#' octahedral strain, \eqn{\bar{\epsilon}_s} (Nádai, 1950) and Lode's parameter
#' for the symmetry of strain \eqn{\nu} (Lode, 1926).
#'
#' @param x accepts the following objects: a two-column matrix where first
#' column is the ratio of maximum strain and
#' intermediate strain (X/Y) and second column is the the ratio of intermediate
#' strain and minimum strain (Y/Z);
#' objects of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"ortensor"` and `"ellipsoid"` objects.
#' Tensor objects can also be lists of such objects (`"ortensor"` and `"ellipsoid"`).
#' @inheritParams Rphi_plot
#' @param labels character. text labels
#' @param add logical. Should data be plotted to an existing plot?
#' @param ... plotting arguments passed to [graphics::points()]
#' @param es.max maximum strain for scaling.
#'
#' @returns a list containing the Lode parameter `lode` and the natural octahedral strain `es`.
#' @family fabric-plot
#' @seealso [lode()] for Lode parameter, and [nadai] for natural octahedral strain.
#' [ellipsoid()] class, [ortensor()] class
#'
#' @details
#' The **amount of strain** related to the natural octahedral unit shear \eqn{\bar{\gamma}_o}
#' is (Nádai, 1963, p.73):
#' \deqn{\bar{\epsilon}_s = \frac{\sqrt{3}}{2} \bar{\gamma}_o}
#' where \eqn{\bar{\gamma}_o} is defined as
#' \deqn{\bar{\gamma}_o = \frac{2}{3} \sqrt{(\bar{\epsilon}_1 - \bar{\epsilon}_2)^2 + (\bar{\epsilon}_2 - \bar{\epsilon}_3)^2 + (\bar{\epsilon}_3 - \bar{\epsilon}_1)^2}}
#' and \eqn{\bar{\epsilon}} is the natural strain (\eqn{\bar{\epsilon} = \log{1+\epsilon}})
#' and \eqn{\epsilon} is the conventional strain given by \eqn{\epsilon = \frac{l-l_0}{l_0}}
#' where \eqn{l} and \eqn{l_0} is the length after and before the strain, respectively (Nádai, 1959, p.70).
#' The amount of strain \eqn{\bar{\epsilon}_s} is directly proportional to the
#' amount of mechanical work applied in the coaxial component of strain.
#'
#' The **symmetry of strain** is defined by Lode’s (1926, p.932) ratio (\eqn{\nu}):
#' \deqn{\nu = \frac{2 \bar{\epsilon}_2 - \bar{\epsilon}_1 - \bar{\epsilon}_3}{\bar{\epsilon}_1 - \bar{\epsilon}_3}}
#' The values range between -1 and +1, where -1 gives constriction, 0 gives
#' plane strain, and +1 gives flattening.
#'
#' @note Hossack (1968) was the first one to use this graphical representation
#' of 3D strain and called the plot "Strain plane plot"
#'
#' @name hsu_plot
#'
#' @references
#' Lode, W. (1926). Versuche über den Einfluß der mittleren Hauptspannung auf
#' das Fließen der Metalle Eisen. Kupfer und Nickel. Zeitschrift Für Physik,
#' 36(11–12), 913–939. \doi{10.1007/BF01400222}
#'
#' Nádai, A. (1950). Theory of flow and fracture of solids. McGraw-Hill.
#'
#' Hsu, T. C. (1966). The characteristics of coaxial and non-coaxial strain
#' paths. Journal of Strain Analysis, 1(3), 216–222. \doi{10.1243/03093247V013216}
#'
#' Hossack, J. R. (1968). Pebble deformation and thrusting in the Bygdin area
#' (Southern Norway). Tectonophysics, 5(4), 315–339. \doi{10.1016/0040-1951(68)90035-8}
#'
#' @examples
#' # default
#' R_XY <- holst[, "R_XY"]
#' R_YZ <- holst[, "R_YZ"]
#' hsu_plot(cbind(R_XY, R_YZ), col = "#B63679", pch = 16, type = "b")
#'
#' # orientation data
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' hsu_plot(x, labels = "x")
#'
#' set.seed(20250411)
#' y <- rvmf(100, mu = mu, k = 20)
#' hsu_plot(ortensor(y), labels = "y", col = "red", add = TRUE)
#'
#' # ellipsoid objects
#' hossack_ell <- lapply(seq.int(nrow(hossack1968)), function(i) {
#'   ellipsoid_from_stretch(hossack1968[i, 3], hossack1968[i, 2], hossack1968[i, 1])
#' })
#' hsu_plot(hossack_ell, col = "#B63679", pch = 16)
NULL

#' @rdname hsu_plot
#' @export
hsu_plot <- function(x, ...) UseMethod("hsu_plot")

#' @rdname hsu_plot
#' @export
hsu_plot.ortensor <- function(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...) {
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

  es <- nadai(x_eigen)
  es.max <- if (is.null(es.max)) max(es) else es.max

  lode <- lode(x_eigen)

  if (isFALSE(add)) {
    hsu_plot.default(cbind(0, 0), es.max = es.max)
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

#' @rdname hsu_plot
#' @export
hsu_plot.spherical <- function(x, ...) {
  ortensor(x) |>
    hsu_plot.ortensor(...)
}

#' @rdname hsu_plot
#' @export
hsu_plot.ellispoid <- function(x, ...) {
  as.ortensor(x) |>
    hsu_plot.ortensor(...)
}


#' @rdname hsu_plot
#' @export
hsu_plot.default <- function(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...) {
  R_XY <- x[, 1]
  R_YZ <- x[, 2]

  R_XZ <- R_XY * R_YZ
  R <- es <- 1 / sqrt(3) * sqrt(log(R_XY)^2 + log(R_YZ)^2 + log(1 / R_XZ)^2) # Nadai, 1963

  K <- log(R_XY) / log(R_YZ) # Hossack 1968
  lode <- (1 - K) / (1 + K)

  # Set plot limits

  rmax <- if (is.null(es.max)) max(R) * 1.1 else es.max
  rseq <- pretty(c(0, rmax))
  rmax2 <- max(rseq)

  if (isFALSE(add)) {
    plot(0, 0,
      type = "n", asp = 1,
      xlim = c(-1, 1) * rmax2 * cos(pi / 2 + pi / 6),
      ylim = c(0, rmax2 * 1.1),
      axes = FALSE, xlab = "", ylab = "", main = main
    )

    # Draw wedge boundaries (±30°)
    graphics::segments(0, 0, rmax2 * cos(pi / 2 + pi / 6), rmax2 * sin(pi / 2 + pi / 6),
      lty = 1, col = "black"
    )
    graphics::segments(0, 0, rmax2 * cos(pi / 2 - pi / 6), rmax2 * sin(pi / 2 - pi / 6),
      lty = 1, col = "black"
    )


    graphics::segments(0, 0, rmax2 * cos(pi / 2 - pi / 6 / 2), rmax2 * sin(pi / 2 - pi / 6 / 2),
      lty = 1, col = "lightgray"
    )

    graphics::segments(0, 0, rmax2 * cos(pi / 2 + pi / 6 / 2), rmax2 * sin(pi / 2 + pi / 6 / 2),
      lty = 1, col = "lightgray"
    )

    # Draw vertical plane strain line
    graphics::segments(0, 0, 0, rmax2, lty = 1, col = "lightgray")


    # Draw cropped concentric arcs (strain magnitude ticks)
    arc_theta <- seq(pi / 2 - pi / 6, pi / 2 + pi / 6, length.out = 200)
    for (rr in rseq) {
      col_grd <- if (rr == rmax2) "black" else "lightgray"
      lwd_grd <- if (rr == rmax2) 1 else 0.5

      graphics::lines(rr * cos(arc_theta), rr * sin(arc_theta), col = col_grd, lwd = lwd_grd)
      # text(rr, 0, labels = format(rr, digits = 2), pos = 4, col = "grey40", cex = 0.8)
    }

    # Lode parameters labels
    for (rr in rseq[-1]) {
      graphics::text(cos(pi / 2 - pi / 6 * 1.05) * rr, sin(pi / 2 - pi / 6 * 1.05) * rr, rr, adj = 1)
    }

    # Strain magnitude labels
    for (vr in seq(-1, 1, by = 0.5)) {
      graphics::text(rmax2 * cos(pi / 2 + pi / 6 * vr) * 1.05, rmax2 * sin(pi / 2 - pi / 6 * vr) * 1.05, vr, adj = 0.5)
    }


    graphics::text(0, rmax2 * 1.1, expression(bold(nu)), adj = 0.5, col = "black", font = 2)
    graphics::text(rmax2 * cos(pi / 2 - pi / 6 * 1.3) / 2, rmax2 * sin(pi / 2 - pi / 6 * 1.3) / 2, expression(bold(bar(epsilon[s]))), adj = .5, col = "black", font = 2)

    # Labels
    graphics::text(0, rmax2 * 0.8, "Plane strain", adj = 0.5, col = "grey70", srt = 90)
    graphics::text(rmax2 * cos(pi / 2 + pi / 6 * 0.9) * 0.8, rmax2 * sin(pi / 2 + pi / 6 * 0.9) * 0.8, "Flattening", adj = 0.5, col = "grey70", srt = 60)
    graphics::text(rmax2 * cos(pi / 2 - pi / 6 * 0.9) * 0.8, rmax2 * sin(pi / 2 - pi / 6 * 0.9) * 0.8, "Constriction", adj = 0.5, col = "grey70", srt = -60)
  }

  # Data points

  # Map Lode parameter (-1..1) to angle in radians (-30°..+30° around vertical)
  theta <- lode * (pi / 6) # -1 -> -30°, 0 -> 0° (vertical), +1 -> +30°

  # Shift so plane strain = vertical (pi/2)
  theta_shift <- theta + pi / 2

  # Cartesian coordinates
  x <- R * cos(theta_shift)
  y <- R * sin(theta_shift)
  graphics::points(x, y, ...)

  invisible(list(lode = lode, es = R))
}

#' @rdname hsu_plot
#' @export
hsu_plot.list <- function(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...) {
  x_eigen <- lapply(x, principal_stretch)
  es <- sapply(x_eigen, nadai)
  es.max <- if (is.null(es.max)) max(es) else es.max
  lode <- sapply(x_eigen, lode)

  if (isFALSE(add)) {
    hsu_plot.default(cbind(0, 0), es.max = es.max, main = main)
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


#' Flinn Diagram
#'
#' Plots the strain ratios X/Y against Y/Z and shows the strain intensity and
#' the strain symmetry after Flinn (1965)
#'
#' @inheritParams hsu_plot
#' @param R.max numeric. Maximum aspect ratio for scaling.
#' @param log logical. Whether the axes should be in logarithmic scale.
#'
#' @returns list. Relative magnitudes of X, Y and Z (Z=1).
#'
#' @family fabric-plot
#' @name flinn_plot
#' @seealso [ellipsoid()] class, [ortensor()] class, [flinn()] for Flinn's
#' strain parameters.
#'
#' @details **Strain symmetry** (Flinn 1965):
#' \deqn{k = \frac{s_1/s_2 - 1}{s_2/s_3 - 1}}
#' where \eqn{s_1 \geq s_2 \geq s_3} the semi-axis lengths of the ellipsoid.
#' The value ranges from 0 to \eqn{\infty}, and is 0 for oblate ellipsoids
#' (flattening), 1 for plane strain and  \eqn{\infty} for prolate ellipsoids (constriction).
#'
#' and **strain intensity** (Flinn 1965):
#' \deqn{d = \sqrt{(s_1/s_2 - 1)^2 + (s_2/s_3 - 1)^2}}
#'
#'
#' @references Flinn, D. (1965). On the Symmetry Principle and the Deformation
#' Ellipsoid. Geological Magazine, 102(1), 36–45. \doi{10.1017/S0016756800053851}
#'
#' @examples
#' data(holst)
#' R_XY <- holst[, "R_XY"]
#' R_YZ <- holst[, "R_YZ"]
#' flinn_plot(cbind(R_XY, R_YZ), log = FALSE, col = "#B63679", pch = 16)
#' flinn_plot(cbind(R_XY, R_YZ), log = TRUE, col = "#B63679", pch = 16, type = "b")
#'
#' # ellipsoid objects
#' hossack_ell <- lapply(seq.int(nrow(hossack1968)), function(i) {
#'   ellipsoid_from_stretch(hossack1968[i, 3], hossack1968[i, 2], hossack1968[i, 1])
#' })
#' flinn_plot(hossack_ell, col = "#B63679", pch = 16, log = TRUE)
#'
#' set.seed(20250411)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 1)
#' flinn_plot(x, R.max = 2)
#'
#' set.seed(20250411)
#' y <- rvmf(100, mu = mu, k = 20)
#' flinn_plot(ortensor(y), col = "red", R.max = 2, add = TRUE)
NULL

#' @rdname flinn_plot
#' @export
flinn_plot <- function(x, main = "Flinn diagram", R.max = NULL, log = FALSE, add = FALSE, ...) UseMethod("flinn_plot")

#' @rdname flinn_plot
#' @export
flinn_plot.default <- function(x, main = "Flinn diagram", R.max = NULL, log = FALSE, add = FALSE, ...) {
  R_XY <- x[, 1]
  R_YZ <- x[, 2]

  if (log) {
    xlab <- "ln(Y/Z)"
    ylab <- "ln(X/Y)"
    R_XY <- log(R_XY)
    R_YZ <- log(R_YZ)
    R.min <- 0
  } else {
    xlab <- "Y/Z"
    ylab <- "X/Y"
    R.min <- 1
  }

  if (is.null(R.max)) {
    R.max <- max(c(R_XY, R_YZ)) * 1.05
  }

  if (isFALSE(add)) {
    plot(R.min, R.min,
      type = "n",
      asp = 1,
      xlim = c(R.min, R.max),
      ylim = c(R.min, R.max),
      # log = log,
      axes = FALSE,
      xaxs = "i", yaxs = "i",
      xlab = xlab,
      ylab = ylab,
      main = main
    )

    graphics::abline(0, b = 1, col = "grey30", lwd = .75)

    rbreaks <- pretty(c(R.min, R.max))


    for (i in c(.2, .5, 2, 5)) {
      graphics::abline(a = 0, b = i, col = "lightgray", lwd = .5, lty = 2)
    }


    if (log) {
      arc_theta <- seq(0, pi / 2, length.out = 200)
      for (rr in rbreaks) {
        graphics::lines(rr * cos(arc_theta), rr * sin(arc_theta), col = "lightgray", lwd = .5, lty = 2)
      }
    } else {
      for (rr in rbreaks) {
        # graphics::abline(a = rr, b = -1, col = "grey80", lwd = .5, lty = 2)
        graphics::lines(c(rr, R.min), c(R.min, rr), col = "lightgray", lwd = .5, lty = 2)
      }
    }

    graphics::axis(side = 1, at = rbreaks, labels = rbreaks)
    graphics::axis(side = 2, at = rbreaks, labels = rbreaks)

    graphics::text(R.max / 2, R.max / 2, "Plane strain", adj = 0.5, col = "grey70", srt = 45)
    graphics::text(R.max * .75, R.max * .25, "Flattening", adj = 0.5, col = "grey70")
    graphics::text(R.max * .25, R.max * .75, "Constriction", adj = 0.5, col = "grey70")
  }

  graphics::points(R_YZ, R_XY, ...)

  invisible(
    list(
      X = R_XY * R_YZ,
      Y = R_YZ,
      Z = 1
    )
  )
}

#' @rdname flinn_plot
#' @export
flinn_plot.ortensor <- function(x, ...) {
  x.stretch <- principal_stretch(x)

  a <- sort(x.stretch, TRUE) |> unname()

  R_xy <- a[1] / a[2]
  R_yz <- a[2] / a[3]

  flinn_plot.default(cbind(R_xy, R_yz), ...)
}

#' @rdname flinn_plot
#' @export
flinn_plot.ellipsoid <- function(x, ...) {
  flinn_plot.ortensor(as.ortensor(x), ...)
}


#' @rdname flinn_plot
#' @export
flinn_plot.spherical <- function(x, ...) {
  flinn_plot.ortensor(ortensor(x), ...)
}

#' @rdname flinn_plot
#' @export
flinn_plot.list <- function(x, ...) {
  a <- lapply(x, principal_stretch) |>
    sapply(sort, decreasing = TRUE) |>
    t()

  R_xy <- a[, 1] / a[, 2]
  R_yz <- a[, 2] / a[, 3]

  flinn_plot.default(cbind(R_xy, R_yz), ...)
}
