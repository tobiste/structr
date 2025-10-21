#' Confidence ellipse
#'
#' Bootstrapped projected mean with percentile confidence region and hypothesis tests
#'
#' @note The inference is based on percentiles of Mahalanobis distance in the tangent
#' space at the mean of the bootstrapped means. The user should check that the
#' bootstrapped means form a tight ellipsoidal cluster, before taking such a
#' region seriously.
#'
#' @inheritParams geodesic_meanvariance_line
#' @param n integer (10000 by default).
#' @param alpha numeric (0.05 by default).
#' @param res integer. The resolution with which to sample the boundary curve of the (1 - `alpha`) * 100% percentile region.
#' @param isotropic logical. If `TRUE`, forces the inverse covariance to the identity matrix and hence the region to be circular.
#'
#' @returns list.
#' \describe{
#' \item{`center`}{Projected mean of `x`}
#' \item{`cov`}{Inverse covariance matrix in the tangent space, which is the identity if `isotropic` is `TRUE`}
#' \item{`rotation`}{Rotation matrix used}
#' \item{`quantiles`}{Quantiles of Mahalanobis norm}
#' \item{`pvalue.ray`}{For rays: p-value for each ray in `x`, i.e. the fraction of `x` that are farther from `center` than the given ray.}
#' \item{`pvalue.line`}{For lines: p-value for each line in `x`, i.e. the fraction of `x` that are farther from `center` than the given line}
#' \item{`pvalue.line.FUN`,`pvalue.ray.FUN`}{The function to calculate the p-value for a given vector}
#' \item{`angles`}{Angles of the semi-axis of the confidence ellipse (in radians if `x` is an `"Vec3"` object, in degrees if otherwise.)}
#' \item{`ellipse`}{Confidence ellipse given as `"Vec3"` object with `res` vectors}
#' }
#' @export
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. https://doi.org/10.1016/j.jsg.2017.01.002
#' @source geologyGeometry (J.R. Davis)
#'
#' @examples
#' set.seed(20250411)
#' ce <- confidence_ellipse(example_lines, n = 1000, res = 10)
#' # print(ce)
#'
#' # Check how many vectors lie outside quantiles:
#' stats::quantile(ce$pvalue, probs = c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))
#'
#' # Hypothesis testing (reject if p-value < alpha):
#' ce$pvalue.FUN((Line(90, 0)))
confidence_ellipse <- function(x, n = 10000L, alpha = 0.05, res = 512L, isotropic = FALSE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Ray(x) | is.Plane(x))

  xl <- vec_list(x)

  if (is.Ray(x)) {
    ce <- rayBootstrapInference(xl, numBoots = n, alpha = alpha, numPoints = res, doIsotropic = isotropic)
    pvalues <- sapply(ce$us, ce$pvalue)
    pvalue.FUN <- ce$pvalue
  } else {
    ce <- lineBootstrapInference(xl, numBoots = n, alpha = alpha, numPoints = res, doIsotropic = isotropic)
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
  ellipse <- if (res > 0) as.Vec3(do.call(rbind, ce$points)) else NULL

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
  )
}


lineBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls,
    length(ls),
    replace = TRUE
  )), simplify = FALSE)
  bootMean <- lineProjectedMean(boots)
  boots <- lapply(boots, function(u) {
    if (dot(u, bootMean) >= 0) {
      u
    } else {
      -u
    }
  })
  inf <- rayMahalanobisPercentiles(boots, center = bootMean, ...)

  inf$pvalueLine <- function(l) {
    if (is.spherical(l)) l <- unclass(Vec3(l))
    if (dot(l, inf$center) < 0) {
      inf$pvalue(-l)
    } else {
      inf$pvalue(l)
    }
  }
  inf
}

rayBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls,
    length(ls),
    replace = TRUE
  )), simplify = FALSE)
  bootMean <- rayProjectedMean(boots)
  rayMahalanobisPercentiles(boots, bootMean, ...)
}


lineProjectedMean <- function(us) {
  eig <- lineMeanScatter(us)
  eig$vectors[, 1]
}


rayProjectedMean <- function(us) rayNormalized(arithmeticMean(us))


lineMeanScatter <- function(us) {
  tMatrices <- lapply(us, function(u) outer(u, u))
  tMatrix <- arithmeticMean(tMatrices)
  eig <- eigen(tMatrix, symmetric = TRUE)
  eig
}

arithmeticMean <- function(xs, weights = NULL) {
  if (is.null(weights)) {
    Reduce("+", xs) / length(xs)
  } else {
    total <- Reduce("+", mapply(function(w, x, SIMPLIFY = TRUE) {
      w * x
    }, weights, xs))
    total / sum(weights)
  }
}


rayMahalanobisPercentiles <- function(us, center, alpha = 0.05, numPoints = 0, doIsotropic = FALSE) {
  perp <- rayOrthogonalUniform(center)
  rot <- rbind(center, perp, cross(center, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  if (doIsotropic) {
    covarInv <- diag(c(1, 1))
  } else {
    covar <- arithmeticMean(lapply(vs, function(v) {
      outer(v, v)
    }))
    covarInv <- solve(covar)
  }
  norms <- sapply(vs, function(v) {
    sqrt(v %*% covarInv %*% v)
  })
  empiricalCDF <- stats::ecdf(norms)
  f <- function(u) {
    if (is.spherical(u)) u <- unclass(Vec3(u))
    v <- rayTangentVectorFromPoint(u, rot)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  qs <- stats::quantile(norms, probs = c(
    0, 0.25, 0.5, 0.75, 0.95,
    0.99, 1
  ), names = FALSE)
  if (alpha != 0.05 && alpha != 0.01) {
    alpha <- 0.05
  }
  if (alpha == 0.01) {
    q <- qs[[6]]
  } else {
    q <- qs[[5]]
  }
  eig <- eigen(covarInv, symmetric = TRUE)
  radii <- q * eig$values^(-1 / 2)
  result <- list(
    us = us, pvalue = f, center = center, covarInv = covarInv,
    rotation = rot, q000 = qs[[1]], q025 = qs[[2]], q050 = qs[[3]],
    q075 = qs[[4]], q095 = qs[[5]], q099 = qs[[6]], q100 = qs[[7]],
    alpha = alpha, angles = radii
  )
  if (numPoints > 0) {
    circle <- lapply(0:numPoints, function(s) {
      c(cos(s *
        2 * pi / numPoints), sin(s * 2 * pi / numPoints))
    })
    vs <- lapply(circle, function(v) {
      as.numeric(eig$vectors %*%
        (radii * v))
    })
    result$points <- lapply(
      vs, rayPointFromTangentVector,
      rot
    )
    result$alpha <- alpha
  }
  result
}

rayTangentVectorFromPoint <- function(q, rotation) {
  v <- rayLog(rotation[1, ], q)
  w <- as.numeric(rotation %*% v)[2:3]
  w
}

rayLog <- function(p, q) {
  normV <- arcCos(dot(p, q))
  w <- rayNormalized(cross(p, q))
  v <- normV * cross(w, p)
  if (dot(v, q) >= 0) {
    v
  } else {
    -v
  }
}

# vnorm
rayNormalized <- function(v) v / sqrt(sum(v^2))

rayOrthogonalUniform <- function(v) rayOrthogonalProjection(v, rayUniform())

rayOrthogonalProjection <- function(pole, v) rayNormalized(v - pole * dot(v, pole) / dot(pole, pole))

rayUniform <- function(n = NULL) {
  if (is.null(n)) {
    cartesianFromHorizontal(c(stats::runif(1, 0, 2 * pi), stats::runif(
      1,
      -1, 1
    )))
  } else {
    replicate(n, rayUniform(), simplify = FALSE)
  }
}


rayPointFromTangentVector <- function(w, rotation) {
  v <- as.numeric(t(rotation) %*% c(0, w))
  q <- rayExp(rotation[1, ], v)
  q
}

rayExp <- function(p, v) {
  normV <- sqrt(dot(v, v))
  if (normV == 0) {
    p
  } else {
    r <- rbind(p, v / normV, cross(p, v / normV))
    as.numeric(t(r) %*% c(cos(normV), sin(normV), 0))
  }
}
