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
#' @param n_iter integer (10000 by default).
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
#' 96, 65–89. \doi{10.1016/j.jsg.2017.01.002}
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @examples
#' set.seed(20250411)
#' ce <- confidence_ellipse(example_lines, n_iter = 1000, res = 10)
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
