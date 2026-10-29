#' Bootstrapped confidence ellipse
#'
#' Bootstrapped projected mean with percentile confidence region and hypothesis tests
#'
#' @note The inference is based on percentiles of Mahalanobis distance in the tangent
#' space at the mean of the bootstrapped means. The user should check that the
#' bootstrapped means form a tight ellipsoidal cluster, before taking such a
#' region seriously.
#'
#' @inheritParams geodesic_meanvariance_line
#' @param n_iter integer. Number of bootstrap iterations (`10000` by default).
#' @param alpha numeric. Level of significance \eqn{\alpha} (`0.05` by default).
#' @param res integer. The resolution with which to sample the boundary curve of the (1 - `alpha`) * 100% percentile region.
#' @param isotropic logical. If `TRUE`, forces the inverse covariance to the identity matrix and hence the region to be circular.
#'
#' @returns object of class `"confidence_region"`, here with additional arguments:
#' \describe{
#' \item{`center`}{Projected mean of `x`}
#' \item{`cov`}{Inverse covariance matrix in the tangent space, which is the identity if `isotropic` is `TRUE`}
#' \item{`rotation`}{Rotation matrix used}
#' \item{`quantiles`}{Quantiles of Mahalanobis norm}
#' \item{`pvalue.ray`}{For rays: p-value for each ray in `x`, i.e. the fraction of `x` that are farther from `center` than the given ray.}
#' \item{`pvalue.line`}{For lines: p-value for each line in `x`, i.e. the fraction of `x` that are farther from `center` than the given line}
#' \item{`pvalue.line.FUN`,`pvalue.ray.FUN`}{The function to calculate the p-value for a given vector}
#' \item{`angles`}{Angles of the semi-axis of the confidence ellipse (in radians if `x` is a `"Vec3"` object, in degrees if otherwise.)}
#' \item{`ellipse`}{Confidence ellipse given as `"Vec3"` object with `res` vectors}
#' }
#' @export
#'
#' @family confidence
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. \doi{10.1016/j.jsg.2017.01.002}
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @examples
#' set.seed(20250411)
#' ce <- confidence_ellipse(Ray(example_lines), n_iter = 1000, res = 10)
#' # print(ce)
#'
#' # Check how many vectors lie outside quantiles:
#' stats::quantile(ce$pvalue, probs = c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))
#'
#' # Hypothesis testing (reject if p-value < alpha):
#' ce$pvalue.FUN((Line(90, 0)))
confidence_ellipse <- function(x, n_iter = 10000L, alpha = 0.05, res = 512L, isotropic = FALSE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Ray(x) | is.Plane(x))

  xl <- vec_list(x)

  if (is.Ray(x)) {
    ce <- rayBootstrapInference(xl, numBoots = n_iter, alpha = alpha, numPoints = res, doIsotropic = isotropic)
    pvalues <- sapply(ce$us, ce$pvalue)
    pvalue.FUN <- ce$pvalue
    ce$center <- as.vector(ce$center)
  } else {
    ce <- lineBootstrapInference(xl, numBoots = n_iter, alpha = alpha, numPoints = res, doIsotropic = isotropic)
    pvalues <- sapply(ce$us, ce$pvalueLine)
    pvalue.FUN <- ce$pvalueLine
  }

  names_res <- names(ce)
  res_quantiles <- names_res[grepl("q", names_res)]

  angles <- if (is.Vec3(x)) {
    ce$angles
  } else {
    rad2deg(ce$angles)
  }
  ellipse <- if (res > 0) list_vec(ce$points) else NULL

  structure(
    list(
      center = as.Vec3(ce$center) |> Spherical(class(x)[1]),
      cov = ce$covarInv,
      rotation = as.Rotation(ce$rotation),
      quantiles = t(do.call(rbind, ce[res_quantiles])),
      pvalue.FUN = ce$pvalue,
      pvalue = pvalues,
      pvalue.FUN = pvalue.FUN,
      angles = angles,
      ellipse = ellipse
    ),
    class = "confidence_region"
  )
}


#' Confidence region for Fisher statistics
#'
#' Circle of angular radius \eqn{\alpha} (degrees) around the mean vector \eqn{\mu}.
#'
#' @param mu Object of class `spherical`
#' @param alpha numeric. Angle in degrees.
#' @param n integer. Resolution of the cone (number of points for the cone line)
#'
#' @returns object of class `"confidence_region"`
#'
#' @family confidence
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(n = 100, Line(120, 55), k = 10)
#' mu <- sph_mean(x)
#' alpha <- sph_confidence_angle(x)
#'
#' r1 <- fisher_confidence(mu, alpha = alpha)
#'
#' plot(x, col = "grey")
#' points(mu, col = "red", pch = 16, cex = 1.5)
#' stereo_lines(r1, col = "red")
fisher_confidence <- function(mu, alpha, n = 512L) {
  azi <- Line(mu)[, "azimuth"]
  secant <- rotate(mu, Line(azi + 90, 0), -alpha)
  angles <- c(seq(0, 360, length.out = n))
  la <- lapply(angles, function(a) {
    rotate(secant, mu, a)
  })
  do.call(rbind, la) |> Vec3()
  
  structure(
    list(
      center = mu,
      ellipse = la
    ),
    class = "confidence_region"
  )
}

#' Confidence region for Bingham statistics
#'
#' Bingham confidence ellipse
#'
#' @param stats Output of [bingham_statistics()]
#' @inheritParams fisher_confidence
#'
#' @family confidence
#' 
#' @export
#'
#' @returns object of class `"confidence_region"`
#'
#' @examples
#' set.seed(20250411)
#' a <- cov(iris[, 1:3])
#' x <- rbingham(100, a, "Line")
#' r <- bingham_statistics(x)
#' r2 <- bingham_confidence(r)
#' plot(x, col = "grey")
#' points(r$mu, col = "red", pch = 16, cex = 1.5)
#' stereo_lines(r2$ellipse, col = "red")
#'
#' r3 <- bingham_statistics(example_planes)
#' stereoplot()
#' points(example_planes, col = "grey", pch = 16)
#' stereo_confidence(bingham_confidence(r3), col = "red")
bingham_confidence <- function(stats, n = 512) {
  g <- names(stats) %in% c("a95", "gamma")
  gamma <- deg2rad(stats[[names(stats)[g]]])
  mu <- stats$mu[1, ] |> Vec3()

  u <- stats$axes[1, ] |> Vec3()
  v <- stats$axes[2, ] |> Vec3()

  theta <- seq(0, 2 * pi, length.out = 512)
  ct <- cos(theta)
  st <- sin(theta)

  g <- gamma
  # polar equation of an ellipse with semi-axes g[1], g[2]
  denom <- sqrt((g[2L] * ct)^2 + (g[1L] * st)^2)
  rho <- numeric(length(theta))
  pos <- denom > 0
  rho[pos] <- g[1L] * g[2L] / denom[pos]

  res <- outer(cos(rho), mu) + sin(rho) * (outer(ct, u) + outer(st, v))

  
  # as.Vec3(res[, 1L, ])
  
  structure(
    list(
      center = mu,
      ellipse = as.Vec3(res[, 1L, ])
    ),
    class = "confidence_region"
  )
}

#' Axial Cross Moments
#'
#' @details Matrix of unit vectors `x`
#' for Eigenvectors `\nu` is \deqn{
#'
#' m_{ij} = \frac{1}{n} \sum_{k=1}^{n}
#' \left( \mathbf{x}_k \cdot \mathbf{\nu}_i \right)^2
#' \left( \mathbf{x}_k \cdot \mathbf{\nu}_j \right)^2 , \qquad i, j = 1, 2, 3
#' }
#'
#' @inheritParams sph_mean
#' @param nu Three Eigenvectors \eqn{\nu} (object of class `Spherical`). See [ot_eigen()] for details
#'
#' @returns 3 x 3 matrix
#' @export
#'
#' @examples
#' x <- example_planes
#' axial_cross_moments(x, ot_eigen(x)$vectors)
axial_cross_moments <- function(x, nu) {
  # m_ij = mean((x . v_i)^2 (x . v_j)^2)
  # x:    n x 3 matrix of unit vectors
  # v: 3 x 3 matrix of orthogonal unit vectors, one per ROW
  X <- Vec3(x)
  vecs <- Vec3(nu)
  P2 <- tcrossprod(X, vecs)^2 # X %*% t(vecs), squared: n x 3
  crossprod(P2) / nrow(P2) # t(P2) %*% P2 / n: 3 x 3
}

#' Confidence ellipse of an Eigenvector of the orientation tensor
#'
#' @details Large-sample (central limit theorem) approximation: the angle of the
#' sample eigenvector \eqn{\nu_i} toward \eqn{\nu_j} is approximately normal with
#' variance \eqn{m_{ij} / (n (\lambda_i - \lambda_j)^2)}, with
#' \eqn{m_{ij} = \text{mean}((x \cdot \nu_i)^2 (x \cdot \nu_j)^2)}.
#' The semi-angles of the region are
#' \deqn{\arcsin(\sqrt(\chi^2_2(\text{level}) \times m_{ij} / (n (\lambda_i - \lambda_j)^2)))},
#' 90&deg; for equal eigenvalues.
#'
#' The moments make it valid for any distribution of axes; for rotationally
#' symmetric data it is the cone of FLE eq. 6.11-6.13. Under the Bingham model
#' \eqn{m_{ij} / (\lambda_i - \lambda_j)^2} tends to \eqn{1 / (2 \Delta_{ij})}, the
#' semi-axes \eqn{ \sqrt{\big(\chi^2_2(\text{level}) / (2 n \Delta_{ij})\big)}}
#' of Bingham (1974, eq. 7.24).
#'
#' @param lambda numeric. Eigenvalues of the orientation tensor, descending. See [ot_eigen()] for details
#' @param nu Corresponding Eigenvectors. Object of class `Spherical`. See [ot_eigen()] for details
#' @param cross_moments 3x3 matrix (see [axial_cross_moments()])
#' @param n integer. Sample size
#' @param level numeric. Confidence level (\eqn{1-\alpha}). Default is `0.95`
#' @param which integer. Index of the Eigenvector of the ellipse is centered on. 
#' One of `1` (the default), `2` or `3`.
#'
#' @family confidence
#'
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @export
#'
#' @returns list
#'
#' @examples
#' x <- example_planes
#' t <- ot_eigen(x)
#' m <- axial_cross_moments(x, ot_eigen(x)$vectors)
#' r4 <- bingham_ellipse(lambda = t$values, nu = t$vectors, cross_moments = m, n = nrow(x), which = 1)
#' print(r4)
#'
#' stereoplot()
#' points(x, col = "grey", pch = 16)
#' stereo_confidence(bingham_confidence(r4), col = "red")
bingham_ellipse <- function(lambda, nu, cross_moments, n, level = 0.95,
                            which = 1L) {
  # Confidence ellipse of an eigenvector of the orientation tensor
  # (large-sample CLT approximation; see the Python docstring).
  #
  # lambdas:          length-3 vector, eigenvalues of the orientation tensor
  # nu:          3 x 3 matrix, corresponding eigenvectors as ROWS
  # cross_moments: 3 x 3 matrix m_ij, see axial_cross_moments()
  # n:             sample size
  # which:         1, 2 or 3 (1-based)
  V <- Vec3(nu)
  stopifnot(nrow(V) == 3L)

  lambdas <- as.numeric(lambda)
  stopifnot(length(lambda) == 3L)

  M <- cross_moments
  stopifnot(nrow(M) == 3 & ncol(M) == 3)

  if (length(which) != 1L || !which %in% 1:3) {
    stop("`which` must be 1, 2 or 3")
  }

  other <- setdiff(1:3, which)
  chi2 <- stats::qchisq(level, df = 2)
  d <- lambda[which] - lambda[other]
  arg <- sqrt(chi2 * M[which, other] / (n * d^2))
  arg[d == 0] <- Inf # equal eigenvalues
  gamma <- ifelse(arg >= 1, 90, asin(pmin(arg, 1)) * 180 / pi)

    list(
      mu    = Vec3(V[which, ]),
      level = level,
      axes  = Vec3(V[other, ]),
      gamma = gamma
    )
}
