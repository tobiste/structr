# Gray regression --------------------------------------------------------------

#
# Implementation of:
#   Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). "On the least-squares
#   fit of small and great circles to spherically projected orientation
#   data." Mathematical Geology 12(3):173-184.
#
# Data model: a set of p unit vectors (direction cosines) a^i = (a1,a2,a3),
# i = 1..p -- poles to planes, lineation directions, etc. A small circle is
# parameterised by its unit-vector center x = (x1,x2,x3) and angular radius
# x4, via
#     a1 x1 + a2 x2 + a3 x3 = cos(x4)                              (eq. 2)
# The angular residual of point i to this circle is
#     R_i = x4 - acos(a1^i x1 + a2^i x2 + a3^i x3)                 (eq. 4)
# and the least-squares fit minimises r = sum(R_i^2) over (x1,x2,x3,x4)
# subject to x1^2+x2^2+x3^2 = 1 (eq. 3). A great circle is the special case
# x4 = pi/2 fixed, leaving only (x1,x2,x3) free.
#
# The paper solves this by a constrained Newton method: at each iterate,
# it builds the gradient q (eq. 6) and full Hessian H (eq. 7) of r, forms
# the Lagrangian for the unit-length constraint, and solves the resulting
# linear system (eq. 11-14) for the Newton step. This file implements that
# iteration directly, plus the initial-value eigenvector estimator
# (eq. 16-17) and the great-vs-small-circle variance-ratio test (eq. 18 and
# the unnumbered V_r equation in the "Statistical significance" section).
#
# Implementation notes / deviations from the paper text (documented so a
# reader checking against the source can see exactly where this departs):
#
#  - H[4,4] (bottom-right entry of the Hessian, eq. 7) is printed in the
#    paper as the bare scalar "2", unlike every other entry of H and q,
#    which are explicit sums over the p points. Differentiating r =
#    sum(R_i^2) twice with respect to x4 (with R_i = x4 - K_i and K_i not
#    depending on x4) gives d^2(R_i^2)/dx4^2 = 2 for every i, so the sum
#    over p points is 2p, not 2. This is almost certainly a dropped
#    summation sign in the original typesetting (every neighbouring cell
#    in the same matrix is a sum), and 2p is what is used here -- verified
#    against a from-scratch re-derivation of the full q and H matrices
#    (all other entries reproduce the paper's eq. 6-7 exactly), and
#    against numerical recovery of known synthetic circles (see the demo
#    at the bottom of this file).
#
#  - The initial-value formula (eq. 16-17) is implemented as a single
#    (possibly non-centered) scatter matrix M = A'A - p*abar*abar', with
#    abar = colMeans(A) for a small-circle fit and abar = 0 for a great-
#    circle fit (exactly as the paper specifies), rather than assembled
#    term-by-term; this is algebraically identical to eq. (16) but avoids
#    looping over p points in R.
#
#  - The paper's worked "concentric" example (Table 2, last row: a small
#    circle fit to one data set and a great circle fit to another,
#    constrained to share a common center) is implemented generically as
#    fit_concentric_circles(), which accepts any mix of fixed-radius and
#    free-radius point sets sharing one center -- rather than being
#    special-cased to exactly one small + one great circle -- since the
#    paper's own description of the method ("residuals R_i are calculated
#    differently, and the summations in q4 and the fourth row of H apply
#    only to the cleavage [free-radius] data") generalises directly.
#
# Style: base R only (no tidyverse dependency), vectorised over
# observations wherever practical.

## Initial-value estimator (eq. 16-17): least-squares plane through the
## poles, used as the starting point for the Newton iteration below.

# abar: 0-vector for a great-circle fit, colMeans(A) for a small-circle fit
.gray_initial_circle_estimate <- function(A, abar) {
  p <- nrow(A)
  M <- crossprod(A) - p * outer(abar, abar) # eq. (16), vectorised
  ev <- eigen(M, symmetric = TRUE)
  x123 <- ev$vectors[, which.min(ev$values)]
  D <- sum(x123 * abar) # eq. (17)
  list(x123 = x123, D = D)
}

## Core Newton-Lagrange iteration (eq. 4-14), shared by the small- and
## great-circle fits.
##
##   A       : p x 3 matrix of unit direction cosines
##   x4_free : TRUE for a small-circle fit (x4 estimated), FALSE for a
##             great circle (x4 fixed at pi/2)
.gray_fit_circle_core <- function(A, x4_free, x_init, tol, max_iter, verbose) {
  p <- nrow(A)
  x123 <- x_init$x123
  x4 <- if (x4_free) x_init$x4 else pi / 2

  r_prev <- Inf
  converged <- FALSE
  iter_used <- 0L

  for (iter in seq_len(max_iter)) {
    iter_used <- iter

    dots <- pmin(pmax(as.numeric(A %*% x123), -1), 1)
    K <- acos(dots) # angular dist. to center
    R <- x4 - K # eq. (4), angular residuals
    r <- sum(R^2)

    sinK <- sin(K)
    # guard points essentially at the center or its antipode (sin K -> 0);
    # these are degenerate for this parameterisation and contribute ~0
    # weight rather than a division blow-up
    tiny <- abs(sinK) < 1e-9
    sinK[tiny] <- sign(sinK[tiny]) * 1e-9
    sinK[sinK == 0] <- 1e-9
    cosK <- cos(K)

    w_q <- R / sinK # eq. (6) weight
    q123 <- 2 * colSums(A * w_q)

    w_H <- (2 / sinK^2) * (1 + R * (cosK / sinK)) # {ak al} weight, eq. (7)
    H33 <- crossprod(A, A * w_H)

    w_H4 <- 2 / sinK
    H34 <- colSums(A * w_H4) # eq. (7), off-diagonal x4 row/col

    if (x4_free) {
      q4 <- 2 * sum(R)
      q <- c(q123, q4)
      H <- rbind(cbind(H33, H34), c(H34, 2 * p)) # 2p, not 2 -- see file header
      Bj <- c(2 * x123, 0) # dg/dx, g = x1^2+x2^2+x3^2-1
    } else {
      q <- q123
      H <- H33
      Bj <- 2 * x123
    }

    Hinv_q <- solve(H, q)
    Hinv_Bt <- solve(H, Bj)
    g <- sum(x123^2) - 1 # should be ~0 if x123 unit
    lambda <- as.numeric((g - sum(Bj * Hinv_q)) / sum(Bj * Hinv_Bt)) # eq. (14)
    dx <- as.numeric(-(Hinv_q + lambda * Hinv_Bt)) # eq. (13)

    if (x4_free) {
      x123 <- x123 + dx[1:3]
      x4 <- x4 + dx[4]
    } else {
      x123 <- x123 + dx
    }
    x123 <- x123 / sqrt(sum(x123^2)) # re-enforce the unit constraint

    if (verbose) cat(sprintf("iter %2d: r = %.8f, |dx| = %.2e\n", iter, r, sqrt(sum(dx^2))))

    if (abs(r - r_prev) < tol * (abs(r_prev) + tol)) {
      converged <- TRUE
      break
    }
    r_prev <- r
  }

  # final residuals at the converged parameters
  dots <- pmin(pmax(as.numeric(A %*% x123), -1), 1)
  K <- acos(dots)
  R <- x4 - K
  r_final <- sum(R^2)

  list(
    x = x123, x4 = x4, r = r_final, residuals = R, iterations = iter_used,
    converged = converged, p = p
  )
}

## Public fitting functions

#' Least-squares small circle through a set of unit vectors (poles)
#'
#' @param A p x 3 matrix (or data.frame) of unit direction cosines
#' @return an object of class `"circle_fit"`: `center` x (unit vector),
#'   angular radius `x4` (radians), residual sum of squares `r`, per-point
#'   angular `residuals`, and Newton-iteration diagnostics
#' @noRd
.gray_fit_small_circle <- function(A, tol = 1e-10, max_iter = 100, verbose = FALSE) {
  A <- as.matrix(A)
  stopifnot(ncol(A) == 3, nrow(A) >= 4) # q = 3 free params + x4 -> need p > 4
  abar <- colMeans(A)
  init <- .gray_initial_circle_estimate(A, abar)
  init$x4 <- acos(pmin(pmax(init$D, -1), 1))
  fit <- .gray_fit_circle_core(A,
    x4_free = TRUE, x_init = init, tol = tol,
    max_iter = max_iter, verbose = verbose
  )
  structure(c(fit, list(type = "small")), class = "gray_circle_fit")
}

#' Least-squares great circle through a set of unit vectors (poles)
#'
#' Equivalent to fit_small_circle() with the angular radius fixed at
#' pi/2 (90 degrees).
#'
#' @inheritParams .gray_fit_small_circle
#' @noRd
.gray_fit_great_circle <- function(A, tol = 1e-10, max_iter = 100, verbose = FALSE) {
  A <- as.matrix(A)
  stopifnot(ncol(A) == 3, nrow(A) >= 3)
  init <- .gray_initial_circle_estimate(A, abar = c(0, 0, 0))
  fit <- .gray_fit_circle_core(A,
    x4_free = FALSE, x_init = init, tol = tol,
    max_iter = max_iter, verbose = verbose
  )
  structure(c(fit, list(type = "great")), class = "gray_circle_fit")
}

#' Least-Squares Small and Great Circle Fit for Spherical Data
#'
#' Finds the best small and great circles using the algorithm by Gray et al. (1980)
#'
#' @inheritParams sph_mean
#' @param out character. Type of regression. `"small"` for a small-circle fit, `"great"` for the great-circle fit, and `"both"` for small- and great-circle estimates.
#' @param tol numeric. Machine precision. Defaults to `getOption("structr.tol")`.
#' @param max_iter integer. Number of Newton-Lagrange optimizations. Default is `100`
#' @param verbose logical. Whether results of the iteration should be printed during the calculation. Default is `FALSE`
#' @param sigma2 numeric. Known error in degrees for a goodness-of-fit test. Default is `10` &deg;
#'
#' @details
#' ## Goodness-of-fit test
#' \eqn{\chi^2} goodness-of-fit check against an independently known error
#' variance \eqn{\sigma^2}
#' \deqn{E(\sum R_i^2) = (p-q) \sigma^2}, \eqn{q=3} (great
#' circle) or \eqn{q=4} (small circle).
#'
#' ## Test whether a small-circle fit is a significant improvement over the
#' corresponding great-circle fit
#'
#' \deqn{V_r = (p-3)  [(r_g - r_c) / r_c] ~ F(1, p-3)} under the paper's
#' normal-residual assumptions.
#'
#' @references Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). On the
#' least-square fit of small and great circles to spherically projected data.
#' Mathematical Geology, Vol. 12, No. 3, 1980.
#'
#' @seealso [regression_smallcircle()], [regression_greatcircle()]
#'
#' @returns list. \describe{
#' \item{`vec`}{`Line` of `Vec3` object, sthe axis of the circle solution}
#' \item{`cone`}{Half-apical angle of the small-circle in degrees. Alway 90 for a great-circle solution.}
#' \item{`r_squared`}{Residual sum of squares}
#' \item{`residuals`}{Per-point residuals}
#' \item{`iterations`, `converged`, `p`}{Newton-iteration diagnostics}
#' \item{`type`}{character. Type of the solution}
#' \item{`statistic`, `df`, `p_value`}{Statistic, degrees of freedom and p-value of the \eqn{\chi^2} goodness-of-fit test}
#' }
#'
#' If `out = 'both`, the a comparison of both solutions is added to the output:
#' \describe{
#' \item{`V_r`}{Test statistic \eqn{V_r}}
#' \item{`df1, df2`}{Degrees of freedom of the small and the great circle data, respectively}
#' \item{`p_value`}{P-value of the comparison test}
#' }
#' @name regression-gray
#'
#' @examples
#' data("gray_example")
#' 
#' # regression:
#' best_clea <- regression_gray(gray_example[1:8, ], out = "both")
#' best_bedd <- regression_gray(gray_example[9:16, ], out = "both")
#' best_all <- regression_gray(gray_example, out = "both")
#'
#' # R-squared values:
#' best_clea$great$r_squared
#' best_bedd$great$r_squared
#' best_all$great$r_squared
#'
#' stereoplot()
#' points(gray_example[1:8, ], col = "blue")
#' points(gray_example[9:16, ], col = "red", pch = "x")
#'
#' # best for cleavage
#' lines(best_clea$small$vec, best_clea$small$cone, col = "lightblue")
#' lines(best_clea$great$vec, 90, lty = 2, col = "blue")
#'
#' # best for bedding
#' lines(best_bedd$small$vec, best_bedd$small$cone, col = "sienna")
#' lines(best_bedd$great$vec, 90, lty = 2, col = "red")
#'
#' # best for all
#' lines(best_all$small$vec, best_bedd$small$cone, col = "gray70")
#' lines(best_all$great$vec, 90, lty = 2, col = "gray60")
#'
#' # Compare small and great-circle solution:
#' regression_gray(gray_example[1:8, ], out = "both")$comparison
#' regression_gray(gray_example[9:16, ], out = "both")$comparison
#' regression_gray(gray_example, out = "both")$comparison
NULL

#' @rdname regression-gray
#' @export
regression_gray <- function(x, out = c("great", "small", "both"), sigma2 = 10, tol = NULL, max_iter = 100, verbose = FALSE) {
  out <- match.arg(out)

  tol <- tol %||% getOption("structr.tol")

  A <- Vec3(x) |> unclass()
  sigma2_rad <- rad2deg(sigma2)

  if (out %in% c("small", "both")) {
    s <- .gray_fit_small_circle(A, tol, max_iter, verbose)
    sgof <- gray_goodness_of_fit(s, sigma2_rad)
    s1 <- s

    s1$x <- as.Vec3(s1$x)
    if (!is.Vec3(x)) {
      s1$x4 <- rad2deg(s$x4)
      s1$x <- Line(s1$x)
    }
    names(s1) <- c("vec", "cone", "r_squared", "residuals", "iterations", "converged", "p", "type")
    s1 <- append(unclass(s1), sgof)
    
  }
  
  if (out %in% c("great", "both")) {
    g <- .gray_fit_great_circle(A, tol, max_iter, verbose)
    ggof <- gray_goodness_of_fit(g, sigma2_rad)
    g1 <- g

    g1$x <- as.Vec3(g$x)
    if (!is.Vec3(x)) {
      g1$x4 <- rad2deg(g$x4)
      g1$x <- Line(g1$x)
    }

    names(g1) <- c("vec", "cone", "r_squared", "residuals", "iterations", "converged", "p", "type")
    g1 <- append(unclass(g1), ggof)
  }

  if (out == "small") {
    return(s1)
  } else if (out == "great") {
    return(g1)
  } else {
    list(
      great = g1, small = s1, comparison = gray_compare_small_vs_great(s, g)
    )
  }
}

#' @rdname regression-gray
#' @export
regression_smallcircle_gray <- function(x, sigma2 = 10, tol = NULL, max_iter = 100, verbose = FALSE){
  tol <- tol %||% getOption("structr.tol")
  
  A <- Vec3(x) |> unclass()
  sigma2_rad <- rad2deg(sigma2)
  
    s <- .gray_fit_small_circle(A, tol, max_iter, verbose)
    sgof <- gray_goodness_of_fit(s, sigma2_rad)

    s$x <- as.Vec3(s$x)
    if (!is.Vec3(x)) {
      s$x4 <- rad2deg(s$x4)
      s$x <- Line(s$x)
    }
    names(s) <- c("vec", "cone", "r_squared", "residuals", "iterations", "converged", "p", "type")
    append(unclass(s), sgof)
}

#' @rdname regression-gray
#' @export
regression_greatcircle_gray <- function(x, sigma2 = 10, tol = NULL, max_iter = 100, verbose = FALSE){
  tol <- tol %||% getOption("structr.tol")
  
  A <- Vec3(x) |> unclass()
  sigma2_rad <- rad2deg(sigma2)
  
  g <- .gray_fit_great_circle(A, tol, max_iter, verbose)
  ggof <- gray_goodness_of_fit(g, sigma2_rad)

  g$x <- as.Vec3(g$x)
  if (!is.Vec3(x)) {
    g$x4 <- rad2deg(g$x4)
    g$x <- Line(g$x)
  }
  
  names(g) <- c("vec", "cone", "r_squared", "residuals", "iterations", "converged", "p", "type")
  append(unclass(g), ggof)
}


## -----------------------------------------------------------------------
## Statistical tests (Sec. "Statistical significance of the fits obtained")
## -----------------------------------------------------------------------

#' Test whether a small-circle fit is a significant improvement over the
#' corresponding great-circle fit (unnumbered V_r equation, p. 178)
#'
#' \deqn{V_r = (p-3) * [(r_g - r_c) / r_c] ~ F(1, p-3)} under the paper's
#' normal-residual assumptions.
#'
#' @param fit_small a fit from .gray_fit_small_circle()
#' @param fit_great a fit from .gray_fit_great_circle(), on the *same* data
#' @return a list with the statistic \eqn{V_r}, degrees of freedom, and p-value
#' @importFrom stats pf
#' @noRd
gray_compare_small_vs_great <- function(fit_small, fit_great) {
  stopifnot(
    fit_small$type == "small", fit_great$type == "great",
    fit_small$p == fit_great$p
  )
  p <- fit_small$p
  Vr <- (p - 3) * ((fit_great$r - fit_small$r) / fit_small$r)
  pval <- stats::pf(Vr, df1 = 1, df2 = p - 3, lower.tail = FALSE)
  list(
    Vr = Vr, df1 = 1, df2 = p - 3, p_value = pval,
    interpretation = if (pval < 0.05) {
      "small circle is a significant improvement over the great circle (p<0.05)"
    } else {
      "no significant improvement over the great circle (p>=0.05)"
    }
  )
}

#' Chi-squared goodness-of-fit check against an independently known error
#' variance sigma^2 (eq. 18): E(sum R_i^2) = (p-q)*sigma^2, q=3 (great
#' circle) or q=4 (small circle).
#' @noRd
#' @importFrom stats pchisq
gray_goodness_of_fit <- function(fit, sigma2) {
  q <- if (fit$type == "small") 4 else 3
  df <- fit$p - q
  stat <- fit$r / sigma2
  list(statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE))
}

## -----------------------------------------------------------------------
## Concentric fit: one shared center, some point sets fit with a fixed
## radius (e.g. a great circle, radius pi/2) and (typically) one point
## set fit with a free radius (a small circle) -- generalises the
## paper's worked bedding+cleavage example (Table 2, last row).
## -----------------------------------------------------------------------

#' Fit several point sets to concentric circles about one shared center
#'
#' @param groups a list of groups, each list(A = matrix, radius = "free")
#'   for a free-radius (small) circle or list(A = matrix, radius = <angle
#'   in radians>) for a fixed-radius circle (e.g. radius = pi/2 for a
#'   great circle). At most one group should have radius = "free" (this
#'   matches the paper's use case and keeps the parameterisation to a
#'   single shared x4).
#' @inheritParams regression_gray2
#' @return a list with the shared center x, the free radius x4 (NA if no
#'   free-radius group was given), total residual sum of squares r
#'   (summed over all groups), and per-group residuals
#' @noRd
fit_concentric_circles <- function(groups, tol = 1e-10, max_iter = 100, verbose = FALSE) {
  free_idx <- which(vapply(groups, function(g) identical(g$radius, "free"), logical(1)))
  stopifnot(length(free_idx) <= 1)
  x4_free <- length(free_idx) == 1

  A_all <- do.call(rbind, lapply(groups, `[[`, "A"))
  abar <- if (x4_free) colMeans(groups[[free_idx]]$A) else c(0, 0, 0)
  init <- .gray_initial_circle_estimate(A_all, abar)
  x123 <- init$x123
  x4 <- if (x4_free) acos(pmin(pmax(sum(x123 * abar), -1), 1)) else NA_real_

  r_prev <- Inf
  iter_used <- 0L
  for (iter in seq_len(max_iter)) {
    iter_used <- iter
    q123 <- c(0, 0, 0)
    q4 <- 0
    H33 <- matrix(0, 3, 3)
    H34 <- c(0, 0, 0)
    H44 <- 0
    r_total <- 0
    group_R <- vector("list", length(groups))

    for (gi in seq_along(groups)) {
      A <- groups[[gi]]$A
      target <- if (identical(groups[[gi]]$radius, "free")) x4 else groups[[gi]]$radius
      dots <- pmin(pmax(as.numeric(A %*% x123), -1), 1)
      K <- acos(dots)
      R <- target - K
      group_R[[gi]] <- R
      r_total <- r_total + sum(R^2)

      sinK <- sin(K)
      tiny <- abs(sinK) < 1e-9
      sinK[tiny] <- sign(sinK[tiny]) * 1e-9
      sinK[sinK == 0] <- 1e-9
      cosK <- cos(K)

      w_q <- R / sinK
      q123 <- q123 + 2 * colSums(A * w_q)
      w_H <- (2 / sinK^2) * (1 + R * (cosK / sinK))
      H33 <- H33 + crossprod(A, A * w_H)

      if (identical(groups[[gi]]$radius, "free")) {
        q4 <- q4 + 2 * sum(R)
        w_H4 <- 2 / sinK
        H34 <- H34 + colSums(A * w_H4)
        H44 <- H44 + 2 * nrow(A) # 2p_group -- see file header re. eq. (7)
      }
    }

    if (x4_free) {
      q <- c(q123, q4)
      H <- rbind(cbind(H33, H34), c(H34, H44))
      Bj <- c(2 * x123, 0)
    } else {
      q <- q123
      H <- H33
      Bj <- 2 * x123
    }

    Hinv_q <- solve(H, q)
    Hinv_Bt <- solve(H, Bj)
    g <- sum(x123^2) - 1
    lambda <- as.numeric((g - sum(Bj * Hinv_q)) / sum(Bj * Hinv_Bt))
    dx <- as.numeric(-(Hinv_q + lambda * Hinv_Bt))

    if (x4_free) {
      x123 <- x123 + dx[1:3]
      x4 <- x4 + dx[4]
    } else {
      x123 <- x123 + dx
    }
    x123 <- x123 / sqrt(sum(x123^2))

    if (verbose) cat(sprintf("iter %2d: r_total = %.8f\n", iter, r_total))
    if (abs(r_total - r_prev) < tol * (abs(r_prev) + tol)) break
    r_prev <- r_total
  }

  list(
    x = x123, x4 = x4, r = r_total, group_residuals = group_R,
    iterations = iter_used, groups = groups
  )
}

# Ramsay regression ------------------------------------------------------------


#' Cone or Plane Best Fit of Conically or Cylindrical Disposed Plane Poles
#'
#' Finding the best fit pole of rotation for a given set of points that are
#' assumed to lie on a mutual small or great circle circle using Ramsay 1967 algorithm
#'
#' @param x matrix, where the rows are the observations and the columns are the coordinates of points
#' @importFrom dplyr mutate summarise
#' @references Ramsay, 1967, p. 18-21
#' @returns numeric vector with
#' \describe{
#' \item{`x`,`y`,`z`}{Cartesian coordinates of best fit pole of plane or cone axis,}
#' \item{`e`}{residual of the sum of square of the deviations of the observed poles to the planes from the best fit pole, and}
#' \item{`K`}{(only for cones) half apical angle of best fit cone (in radians).}
#' }
#' @name best_pole
#'
#' @references Ramsay, J. G. (1967). Folding and Fracturing of Rocks. McGraw-Hill.
#'
#' @examples
#' \dontrun{
#' # example from Ramsay, 1967, p. 20
#' x <- rbind(
#'   c(-67, -31, -71),
#'   c(-62, -53, -50),
#'   c(-62, -75, -34),
#'   c(-58, 85, -34),
#'   c(-79, 40, -52),
#'   c(90, 14, -75),
#'   c(80, 10, 90)
#' ) |> acoscartesian_to_cartesian()
#' regression_cone_ramsay(x) # expect: c(0.856, -0.157, -0.492, NA, 1.56207)
#' regression_plane_ramsay(x) # expect: c(0.852, -0.154, -0.502, 1-1.002)
#' }
NULL

regression_cone_ramsay <- function(x) {
  # ensure x is a matrix
  x <- as.matrix(x)
  l <- x[, 1]
  m <- x[, 2]
  n <- x[, 3]

  # precompute sums directly (avoids mutate/summarise)
  l2 <- sum(l^2)
  m2 <- sum(m^2)
  lm <- sum(l * m)
  ln <- sum(l * n)
  mn <- sum(m * n)
  l_sum <- sum(l)
  m_sum <- sum(m)
  n_sum <- sum(n)
  N <- nrow(x)

  # construct matrices directly
  D <- matrix(c(
    l2, lm, l_sum,
    lm, m2, m_sum,
    l_sum, m_sum, N
  ), nrow = 3, byrow = TRUE)

  Da <- matrix(c(
    -ln, lm, l_sum,
    -mn, m2, m_sum,
    -n_sum, m_sum, N
  ), nrow = 3, byrow = TRUE)

  Db <- matrix(c(
    l2, -ln, l_sum,
    lm, -mn, m_sum,
    l_sum, -n_sum, N
  ), nrow = 3, byrow = TRUE)

  Dc <- matrix(c(
    l2, lm, -ln,
    lm, m2, -mn,
    l_sum, m_sum, -n_sum
  ), nrow = 3, byrow = TRUE)

  # determinants
  detD <- det(D)
  A <- det(Da) / detD
  B <- det(Db) / detD
  C <- det(Dc) / detD

  # direction cosines
  cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta <- B * cos_gamma
  cos_K <- -C * cos_gamma

  # half apical angle
  K <- acos(cos_K)

  cart <- c(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
  e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2

  c(cart, e = 1 - e, K = K)
}


#' @rdname best_pole
regression_cone_ramsay2 <- function(x) {
  # ensure matrix
  x <- as.matrix(x)
  l <- x[, 1]
  m <- x[, 2]
  n <- x[, 3]

  # compute sums directly
  l2 <- sum(l^2)
  m2 <- sum(m^2)
  lm <- sum(l * m)
  ln <- sum(l * n)
  mn <- sum(m * n)

  # precompute
  t <- 1 / (l2 * m2 - lm^2)
  A <- (lm * mn - ln * m2) * t
  B <- (lm * ln - mn * l2) * t

  # direction cosines
  cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
  cos_alpha <- A * cos_gamma
  cos_beta <- B * cos_gamma

  # Cartesian coordinates
  cart <- c(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
  e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2

  c(cart, e = 1 - e)
}

# Geodesic regression ----------------------------------------------------------


#' Least-Squares Geodesic Regression of Small and Great Circles for Spherical Data
#'
#' Fits a great circle arc to a set of lines, based on an independent scalar variable.
#'
#' @inheritParams sph_mean
#' @param val A vector of real numbers. Values of independent scalar variables. Same length as `x`.
#' @param iterations A real number (positive integer). A bound on the number of numerical optimization iterations allowed.
#' @param n_points A real number (positive integer). The number of points along the regression curve requested.
#' @param num_seeds A real number (positive integer)
#'
#' @returns list. \describe{
#' \item{`vec`}{axis of best-fit great- or small-circle. Class as `x`}
#' \item{`cone`}{half-apical angle of the best-fit small-circle (in radians if `x` is of class `"Vec3"`, degrees otherwise)}
#' \item{`convergence`}{An integer code. `0` indicates successful completion.}
#' \item{`min_eigenvalue`}{Positive if successful completion}
#' \item{`r_squared`}{\eqn{R^2} value of the regression}
#' \item{`points`}{A set of Vec3 points for the regression line}
#' \item{`range`}{Range of the regression line (in radians if `x` is of class `"Vec3"`, degrees otherwise)}
#' }
#' @name best_fit
#' @family geodesic-regression
#'
#' @references Davis, J.R. and Titus, S.J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology, 96,
#' 65–89. \doi{10.1016/j.jsg.2017.01.002}
#'
#' @examples
#' set.seed(20250411)
#' data("gray_example")
#' bestgc_clea <- regression_greatcircle(gray_example[1:8, ])
#' bestgc_bedd <- regression_greatcircle(gray_example[9:16, ])
#'
#' bestsc_clea <- regression_smallcircle(gray_example[1:8, ])
#' bestsc_bedd <- regression_smallcircle(gray_example[9:16, ])
#' bestsc_all <- regression_smallcircle(gray_example)
#'
#' stereoplot()
#' points(gray_example[1:8, ], col = "blue")
#' points(gray_example[9:16, ], col = "red", pch = "x")
#'
#' # best for cleavage
#' lines(bestsc_clea$vec, bestsc_clea$cone, col = "lightblue")
#' points(bestsc_clea$vec, col = "lightblue", pch = 16)
#' lines(bestgc_clea$vec, lty = 2, col = "blue")
#' points(bestgc_clea$vec, col = "blue", pch = 17)
#'
#' # best for bedding
#' lines(bestsc_bedd$vec, bestsc_bedd$cone, col = "sienna")
#' points(bestsc_bedd$vec, col = "sienna", pch = 16)
#' lines(bestgc_bedd$vec, lty = 2, col = "red")
#' points(bestgc_bedd$vec, col = "red", pch = 17)
NULL

#' @rdname best_fit
#' @export
regression_greatcircle <- function(x, val = seq_len(nrow(x)), iterations = 1000L, n_points = 100L) {
  ls <- vec_list(x)
  n <- length(ls)
  val_num <- as.numeric(val)

  # Let l0 be the l whose x is closest to zero.
  l0 <- ls[[which.min(val_num^2)]]

  # Define the function to be minimized.
  e <- function(wb) {
    a <- rotExp(rotAntisymmetricFromVector(wb[1:3]))
    w <- wb[[4L]]

    # Vectorised trig: compute all unit vectors at once as a 3×n matrix
    cos_w <- cos(w * val_num)
    sin_w <- sin(w * val_num)
    u_mat <- rbind(cos_w, sin_w, 0)
    au_mat <- a %*% u_mat
    # f <- function(i) {
    #   u <- c(cos(wb[[4]] * val[[i]]), sin(wb[[4]] * val[[i]]), 0)
    #   lineDistance(ls[[i]], as.numeric(a %*% u))^2
    # }
    # d <- vapply(1:n, f, numeric(1))
    d <- vapply(
      seq_len(n), function(i) lineDistance(ls[[i]], au_mat[, i]),
      numeric(1L)
    )

    sum(d^2) / (2 * n)
  }

  # Find the minimum, using the constant geodesic l0 as the seed.
  seed <- c(0, 0, 0, 0)
  solution <- stats::optim(seed, e, hessian = TRUE, control = list(maxit = iterations))

  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric = TRUE, only.values = TRUE)$values
  rot <- rotExp(rotAntisymmetricFromVector(solution$par[1:3]))
  a <- solution$par[[4]]

  rSq <- if (is.Line(x) | is.Plane(x)) {
    1 - solution$value / lineVariance(ls, lineProjectedMean(ls))
  } else if (is.Ray(x) | is.Vec3(x)) {
    1 - solution$value / rayVariance(ls, rayProjectedMean(ls))
  }

  rot_vec <- as.Vec3(rot[, 3])

  result <- list(
    vec = Spherical(rot_vec, class(x)[1L]),
    range = a,
    # rotation = rot,
    convergence = solution$convergence,
    min_eigenvalue = min(eigvals),
    r_squared = rSq
  )

  if (!is.Vec3(x)) result$range <- rad2deg(result$range)

  result$prediction <- function(x) {
    pts <- as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
    as.Vec3(pts)
  }

  if (n_points >= 1) {
    result$points <- lapply(
      seq(from = min(val), to = max(val), length.out = (n_points + 1)),
      result$prediction
    ) |>
      list_vec()
  }

  return(result)
}

#' @rdname best_fit
#' @export
regression_smallcircle <- function(x, val = seq_len(nrow(x)), num_seeds = 5L, iterations = 1000L, n_points = 100L) {
  xv <- vec_list(x)

  if (is.Ray(x) | is.Vec3(x)) {
    res <- rayRegressionSmallCircleRescaled(val, xv, numSeeds = num_seeds, numSteps = iterations, numPoints = n_points)
  } else {
    res <- lineRegressionSmallCircleRescaled(val, xv, numSeeds = num_seeds, numSteps = iterations, numPoints = n_points)
  }

  vec <- as.Vec3(res$pole)
  points <- list_vec(res$points)
  cone <- angle(points, vec)[1L]

  if (!is.Vec3(x)) cone <- rad2deg(cone)
  if (!is.Vec3(x)) res$angle <- rad2deg(res$angle)

  list(
    vec = Spherical(vec, class(x)[1L]),
    cone = cone,
    convergence = res$error,
    min_eigenvalue = res$minEigenvalue,
    r_squared = res$rSquared,
    points = points,
    range = res$angle
  )
}


## Permutation

#' Permutation Tests of Regressions
#'
#' `perm_rsq` retrieves the \eqn{R^2} values from successful `n_perm` permutations
#' of a regression function.
#' `perm_rsq_pvalue` post-processes this vector to count how many of its entries
#' exceed the original \eqn{R^2} value. That fraction is a p-value for the null
#' hypothesis that the regression result is meaningless.
#'
#' @param n_perm A real number (positive integer). The number of permutations to try.
#' @param FUN An R function. Must return a result list with fields `convergence`,
#' `min_eigenvalue`, `r_squared` Typically [regression_greatcircle()] or a
#' similar regression function.
#' @param x A real vector. The values of the independent variable. Assumed to be
#' the first argument passed to `FUN`.
#' @param ... Other arguments passed to `FUN`, after `x`.
#' @param r_perm vector of permutated \eqn{R^2} values
#' @param r original \eqn{R^2} value
#'
#' @return `perm_rsq` returns a a real vector. The maximum length is `n_perms`.
#' Often the length is less than `n_perms`, because the regression failed (as
#' signaled by error != 0 or `min_eigenvalue` <= 0).
#'
#' `perm_rsq_pvalue` returns the p-value.
#'
#' @name perm-rsq
#' @family geodesic-regression
#' @importFrom future.apply future_vapply
#'
#' @examples
#' set.seed(20250411)
#' data("gray_example")
#'
#' # original regression
#' bestgc_clea <- regression_greatcircle(gray_example[1:8, ])
#'
#' # permutation test
#' pr <- perm_rsq(100, FUN = regression_greatcircle, x = gray_example[1:8, ])
#'
#' # p-value
#' perm_rsq_pvalue(pr, bestgc_clea)
NULL

#' @name perm-rsq
#' @export
perm_rsq <- function(n_perm, FUN, x, ...) {
  # print("Percentage of permutations complete:")
  g <- function(i) {
    # print(paste0(i / n_perm * 100, "%"))
    newXs <- sample_spherical(x, size = nrow(x))
    regr <- FUN(x = newXs)
    c(regr$convergence, regr$min_eigenvalue, regr$r_squared)
  }

  # perms <- vapply(seq_len(n_perm), g, numeric(3))
  perms <- future.apply::future_vapply(
    seq_len(n_perm), g,
    FUN.VALUE = numeric(3),
    future.seed = TRUE
  )

  success <- perms[1L, ] == 0L & perms[2L, ] > 0
  perms[3L, success]
}

#' @name perm-rsq
#' @export
perm_rsq_pvalue <- function(r_perm, r) {
  if (is.list(r)) r <- r$r_squared
  sum(r_perm > r) / length(r_perm)
}
