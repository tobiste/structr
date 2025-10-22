#' Least-square fit of small and great circles to spherically projected data
#'
#' Finds the best small and great circles using the algorithm by Gray et al. (1980)
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`.
#'
#' @references Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). On the
#' least-square fit of small and great circles to spherically projected data.
#' Mathematical Geology, Vol. 12, No. 3, 1980.
#' @noRd
#' 
#' @seealso [regression_smallcircle()], [regression_greatcircle()]
#'
#' @returns list. `axis_c` is the axis of the small-circle cone, `axis_g` is the axis of the great-circle,
#' `cone_angle` is the half-apical angle of the cone, `r_*` is the residual (in the angle between pole and great/small circle), `E_*` are the errors, and `Vr` is the variance ratio.
#'
#' @examples
#' data("gray_example")
#' best_clea <- regression_gray(gray_example[1:8, ])
#' best_bedd <- regression_gray(gray_example[9:16, ])
#' best_all <- regression_gray(gray_example)
#'
#' stereoplot()
#' points(gray_example[1:8, ], col = "blue")
#' points(gray_example[9:16, ], col = "red", pch = "x")
#'
#' # best for cleavage
#' lines(best_clea$axis_c, best_clea$cone_angle, col = "lightblue")
#' lines(best_clea$axis_g, 90, lty = 2, col = "blue")
#'
#' # best for bedding
#' lines(best_bedd$axis_c, best_bedd$cone_angle, col = "sienna")
#' lines(best_bedd$axis_g, 90, lty = 2, col = "red")
#'
#' # best for all
#' lines(best_all$axis_c, best_bedd$cone_angle, col = "gray70")
#' lines(best_all$axis_g, 90, lty = 2, col = "gray60")
regression_gray <- function(x) {
  # convert from spherical to cartesian coordinates:
  xv <- Vec3(x) #|> unclass()

  p <- nrow(xv)
  g_res <- gray_algorithm(xv, sm = FALSE)
  c_res <- gray_algorithm(xv, sm = TRUE)

  # angle between eigenvector and mean plane/cone
  K_g <- pi / 2 # == acos(g_res$cos_K)
  K_c <- acos(c_res$cos_K)


  # Angle between eigenvector (pole of plane/cone) and data:
  a_g <- angle(g_res$eig3, xv)
  a_c <- angle(c_res$eig3, xv)
  # R <- K - acos(A_eigen3[1] * t(a1) + A_eigen3[2] * t(a2) + A_eigen3[3] * t(a3))
  # R_s <- K_s - acos(c_res$eig3 %*% t(xv))

  # residuals
  R_g <- K_g - a_g
  R_c <- K_c - a_c

  # sum of squares of residuals
  r_g <- sum(R_g^2)
  r_c <- sum(R_c^2)

  E_c <- (p - 2) * r_c
  E_g <- (p - 3) * r_g

  Vr <- (p - 3) * ((r_g - r_c) / r_c)

  axis_c <- as.Vec3(c_res$eig3)
  axis_g <- as.Vec3(g_res$eig3)

  if (is.Line(x) | is.Ray(x) | is.Plane(x)) {
    axis_s <- Line(axis_c)
    axis_g <- Line(axis_g)

    K_c <- rad2deg(K_c)
  }

  return(
    list(
      axis_c = axis_c,
      axis_g = axis_g,
      cone_angle = K_c,
      r_c = r_c,
      r_g = r_g,
      error_c = E_c,
      error_g = E_g,
      variance_ratio = Vr
    )
  )
}


gray_algorithm <- function(x, sm = TRUE) {
  a1 <- x[, 1]
  a2 <- x[, 2]
  a3 <- x[, 3]

  p <- length(a1)
  if (sm) {
    a1_bar <- sum(a1 / p)
    a2_bar <- sum(a2 / p)
    a3_bar <- sum(a3 / p)
  } else {
    a1_bar <- a2_bar <- a3_bar <- 0
  }

  a1_d <- a1 - a1_bar
  a2_d <- a2 - a2_bar
  a3_d <- a3 - a3_bar

  A1 <- rbind(
    sum(a1 * a1_d),
    sum(a2 * a1_d),
    sum(a3 * a1_d)
  )
  A2 <- rbind(
    sum(a1 * a2_d),
    sum(a2 * a2_d),
    sum(a3 * a2_d)
  )
  A3 <- rbind(
    sum(a1 * a3_d),
    sum(a2 * a3_d),
    sum(a3 * a3_d)
  )
  A <- cbind(A1, A2, A3)

  eig3 <- eigen(A, symmetric = TRUE)$vectors[, 3]

  # dot product between eigenvector and mean vector
  cos_K <- eig3[1] * a1_bar + eig3[2] * a2_bar + eig3[3] * a3_bar
  # cos_K <- eig3 %*% t(cbind(a1_bar, a2_bar, a3_bar))
  # cos_K <- vdot(eig3, Vec3(a1_bar, a2_bar, a3_bar))

  return(list(eig3 = as.Vec3(eig3), cos_K = cos_K))
}





#' The cone or plane best fit of conically or cylindrical disposed s-plane poles
#'
#' Finding the best fit pole of rotation for a given set of points that are
#' assumed to lie on a mutual small or great circle circle
#'
#' @param x matrix. Cartesian coordinates of points
#' @importFrom dplyr mutate summarise
#' @references Ramsay, 1967, p. 18-21
#' @returns numeric vector with
#' \describe{
#' \item{`x`,`y`,`z`}{Cartesian coordinates of best fit pole of plane or cone axis,}
#' \item{`e`}{residual of the sum of square of the deviations of the observed poles to the planes from the best fit pole, and}
#' \item{`K`}{(only for cones) half apical angle of best fit cone (in radians).}
#' }
#' @name best_pole
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

#' @rdname best_pole
# best_cone_ramsay <- function(x) {
#   l <- m <- n <- l2 <- m2 <- lm <- ln <- mn <- numeric()
#   xsum <- data.frame(l = x[, 1], m = x[, 2], n = x[, 3]) |>
#     dplyr::mutate(
#       l2 = l^2,
#       m2 = m^2,
#       lm = l * m,
#       ln = l * n,
#       mn = m * n
#     ) |>
#     dplyr::summarise(
#       l = sum(l),
#       m = sum(m),
#       n = sum(n),
#       l2 = sum(l2),
#       m2 = sum(m2),
#       lm = sum(lm),
#       ln = sum(ln),
#       mn = sum(mn)
#     )
#   N <- nrow(xsum)
#
#   D <- Da <- Db <- Dc <- matrix(nrow = 3, ncol = 3)
#   D[1, 1] <- xsum$l2
#   D[1, 2] <- D[2, 1] <- xsum$lm
#   D[1, 3] <- D[3, 1] <- xsum$l
#   D[2, 2] <- xsum$m2
#   D[2, 3] <- D[3, 2] <- xsum$m
#   D[3, 3] <- N
#
#   Da[1, 1] <- -xsum$ln
#   Da[1, 2] <- D[1, 2]
#   Da[1, 3] <- D[1, 3]
#   Da[2, 1] <- -xsum$mn
#   Da[2, 2] <- D[2, 2]
#   Da[2, 3] <- Da[3, 2] <- D[2, 3]
#   Da[3, 1] <- -xsum$n
#   Da[3, 3] <- N
#
#   Db[1, 1] <- D[1, 1]
#   Db[1, 2] <- Da[1, 1]
#   Db[1, 3] <- D[1, 3]
#   Db[2, 1] <- D[1, 2]
#   Db[2, 2] <- Da[2, 1]
#   Db[2, 3] <- D[2, 3]
#   Db[3, 1] <- D[1, 3]
#   Db[3, 2] <- Da[3, 1]
#   Db[3, 3] <- N
#
#   Dc[1, 1] <- D[1, 1]
#   Dc[1, 2] <- D[1, 2]
#   Dc[1, 3] <- Da[1, 1]
#   Dc[2, 1] <- D[1, 2]
#   Dc[2, 2] <- D[2, 2]
#   Dc[2, 3] <- Da[2, 1]
#   Dc[3, 1] <- D[1, 3]
#   Dc[3, 2] <- D[2, 3]
#   Dc[3, 3] <- Da[3, 1]
#
#   A <- det(Da) / det(D)
#   B <- det(Db) / det(D)
#   C <- det(Dc) / det(D)
#
#   # gamma <- -acos((1 + A^2 + B^2)^(-1 / 2))
#   # alpha <- pi - acos(A * (1 + A^2 + B^2)^(-1 / 2))
#   # beta <- -acos(B * (1 + A^2 + B^2)^(-1 / 2))
#   cos_gamma <- 1 / sqrt(1 + A^2 + B^2)
#   cos_alpha <- A * cos_gamma
#   cos_beta <- B * cos_gamma
#   cos_K <- -C * cos_gamma
#
#   # half apical angle
#   K <- acos(cos_K)
#
#   cart <- cbind(x = -cos_alpha, y = -cos_beta, z = -cos_gamma)
#   e <- cos_alpha^2 + cos_beta^2 + cos_gamma^2
#
#   # alpha <- acos(cos_alpha)
#   # beta <- acos(cos_beta)
#   # gamma <- acos(cos_gamma)
#
#   # correct for lower hemisphere and convert to Cartesian coordinates
#   # cart <- cbind(pi-alpha, -beta, -gamma) |>
#   #   tectonicr::rad2deg() |>
#   #   acoscartesian_to_cartesian()
#   # #names(cart) <- NULL
#   # e <- cos(alpha)^2 + cos(beta)^2 + cos(gamma)^2
#
#
#   return(c(cart[, 1], cart[, 2], cart[, 3], "e" = 1 - e, "K" = K))
# }
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


#' Least-square fit of small and great circles to spherically projected data
#'
#' Fits a great circle arc to a set of lines, based on an independent scalar variable.
#' 
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`.
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
#' 
#' @examples
#' data("gray_example")
#' bestgc_clea <- regression_greatcircle(gray_example[1:8, ])
#' bestgc_bedd <- regression_greatcircle(gray_example[9:16, ])
#' bestgc_all <- regression_greatcircle(gray_example)
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
#'
#' # best for all
#' lines(bestsc_all$vec, bestsc_all$cone, col = "gray80")
#' points(bestsc_all$vec, col = "gray80", pch = 16)
#' lines(bestgc_all$vec, lty = 2, col = "gray50")
#' points(bestgc_all$vec, col = "gray50", pch = 17)
NULL

#' @rdname best_fit
#' @export
regression_greatcircle <- function(x, val = seq_len(nrow(x)), iterations = 1000L, n_points = 100L) {
  ls <- vec_list(x)

  # Let l0 be the l whose x is closest to zero.
  l0 <- ls[[which.min(as.numeric(val)^2)]]
  # Define the function to be minimized.
  n <- length(ls)
  e <- function(wb) {
    a <- rotExp(rotAntisymmetricFromVector(wb[1:3]))
    f <- function(i) {
      u <- c(cos(wb[[4]] * val[[i]]), sin(wb[[4]] * val[[i]]), 0)
      lineDistance(ls[[i]], as.numeric(a %*% u))^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic l0 as the seed.
  seed <- c(0, 0, 0, 0)
  solution <- stats::optim(seed, e, hessian = TRUE, control = list(maxit = iterations))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric = TRUE, only.values = TRUE)$values
  rot <- rotExp(rotAntisymmetricFromVector(solution$par[1:3]))
  a <- solution$par[[4]]

  if (is.Line(x) | is.Plane(x)) {
    rSq <- 1 - solution$value / lineVariance(ls, lineProjectedMean(ls))
  } else if (is.Ray(x)) {
    rSq <- 1 - solution$value / rayVariance(ls, rayProjectedMean(ls))
  }

  rot_vec <- as.Vec3(rot[, 3])


  result <- list(
    vec = Spherical(rot_vec, class(x)[1]),
    range = a,
    # rotation = rot,
    convergence = solution$convergence,
    min_eigenvalue = min(eigvals),
    r_squared = rSq
  )

  if (!is.Vec3(x)) result$range <- rad2deg(result$range)


  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }

  if (n_points >= 1) {
    result$points <- do.call(rbind, lapply(
      seq(from = min(val), to = max(val), length.out = (n_points + 1)),
      result$prediction
    ))
  }

  return(result)
}

#' @rdname best_fit
#' @export
regression_smallcircle <- function(x, val = seq_len(nrow(x)), num_seeds = 5L, iterations = 1000L, n_points = 100L) {
  xv <- vec_list(x)
 
  if(is.Ray(x)){
    res <- rayRegressionSmallCircleRescaled(val, xv, numSeeds = num_seeds, numSteps = iterations, numPoints = n_points)
  } else {
    res <- lineRegressionSmallCircleRescaled(val, xv, numSeeds = num_seeds, numSteps = iterations, numPoints = n_points)
  }

  vec <- as.Vec3(res$pole)
  points <- as.Vec3(do.call(rbind, res$points))
  cone <- angle(points,vec)[1]
  
  if (!is.Vec3(x)) cone <- rad2deg(cone)
  if (!is.Vec3(x)) res$angle <- rad2deg(res$angle)
  
  list(
    vec = Spherical(vec, class(x)[1]),
    cone = cone,
    convergence = res$error,
    min_eigenvalue = res$minEigenvalue,
    r_squared = res$rSquared,
    points = points,
    range = res$angle
  )
}

rotMatrixFromAxisAngle <- function(ua) {
  m <- rotAntisymmetricFromVector(ua[1:3])
  diag(c(1, 1, 1)) + sin(ua[[4]]) * m + (1 - cos(ua[[4]])) * (m %*% m)
}

rotAntisymmetricFromVector <- function(v) {
  matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}

lineUniform <- function(n = NULL) {
  if (is.null(n)) {
    lower(rayUniform())
  } else {
    lapply(rayUniform(n), lower)
  }
}

rayRegressionSmallCircle <- function(xs, us, numSeeds = 5, numSteps = 1000, numPoints = 0, angleBound = Inf) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      r <- rotMatrixFromAxisAngle(c(pole, x * phiThetaAlphaTauSigma[[3]]))
      as.numeric(r %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(rayDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value = (2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(
      sphericalFromCartesian(rayUniform())[2:3],
      runif(1, -pi, pi),
      sphericalFromCartesian(rayUniform())[2:3]
    )
    if (angleBound == Inf) {
      solution <- optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B"
      )
    } else {
      solution <- optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B",
        lower = c(-Inf, -Inf, -angleBound, -Inf, -Inf),
        upper = c(Inf, Inf, angleBound, Inf, Inf)
      )
    }
    if (solution$value <= best$value) {
      best <- solution
    }
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric = TRUE, only.values = TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- rayVariance(us, rayProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {
    as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  }
  results <- list(
    error = best$convergence, minEigenvalue = min(eigvals), pole = pole, angle = angle,
    rSquared = rSquared, prediction = pred
  )
  if (numPoints >= 1) {
    ys <- seq(from = min(xs), to = max(xs), length.out = numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

# Warning: Not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
rayRegressionSmallCircleRescaled <- function(xs, ls, angleBound = Inf, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- rayRegressionSmallCircle(
    scales(xs), ls,
    angleBound = (angleBound * (x1 - x0)), ...
  )
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {
    results$rescaledPrediction((x - x0) / (x1 - x0))
  }
  results$angle <- results$angle / (x1 - x0)
  results
}

lineRegressionSmallCircle <- function(xs, us, numSeeds = 5, numSteps = 1000, numPoints = 0, angleBound = Inf) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      angle <- x * phiThetaAlphaTauSigma[[3]]
      as.numeric(rotMatrixFromAxisAngle(c(pole, angle)) %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(lineDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value = (2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(
      sphericalFromCartesian(lineUniform())[2:3],
      runif(1, -pi, pi),
      sphericalFromCartesian(lineUniform())[2:3]
    )
    if (angleBound == Inf) {
      solution <- optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B"
      )
    } else {
      solution <- optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B",
        lower = c(-Inf, -Inf, -angleBound, -Inf, -Inf),
        upper = c(Inf, Inf, angleBound, Inf, Inf)
      )
    }
    if (solution$value <= best$value) {
      best <- solution
    }
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric = TRUE, only.values = TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- lineVariance(us, lineProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {
    as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  }
  results <- list(
    error = best$convergence, minEigenvalue = min(eigvals), pole = pole, angle = angle,
    rSquared = rSquared, prediction = pred
  )
  if (numPoints >= 1) {
    ys <- seq(from = min(xs), to = max(xs), length.out = numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

# Not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
lineRegressionSmallCircleRescaled <- function(xs, us, angleBound = Inf, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- lineRegressionSmallCircle(
    xs = scales(xs), us = us, angleBound = (angleBound * (x1 - x0)), ...
  )
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {
    results$rescaledPrediction((x - x0) / (x1 - x0))
  }
  results$angle <- results$angle / (x1 - x0)
  results
}

# The number interval [a, b] is mapped to the interval [0, 1]. Inputs outside [a, b] are clamped to [a, b].
scale <- function(y, a = 0, b = 1) {
  (y - a) / (b - a)
}

# The number interval range(ys) is mapped to the interval [0, 1].
scales <- function(ys) {
  sapply(ys, scale, min(ys), max(ys))
}

