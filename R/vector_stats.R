#' Mean resultant of a set of vectors
#'
#' @param x numeric. Can be three element vector or a three column array
#' @param w numerical vector of weights the same length as `x` giving the
#' weights to use for elements of `x`.
#' @param mean logical. Whether the mean resultant (`TRUE`) or resultant
#' (`FALSE`, the default) is returned.
#' @returns if `mean==TRUE`, mean resultant is returned
#' and numeric otherwise.
#' @noRd
#' @examples
#' x <- rvmf(100, mu = Line(120, 50), k = 5) |> Vec3() |> unclass()
#' vresultant(x, mean = FALSE)
#' vresultant(x, mean = TRUE)
vresultant <- function(x, w = NULL, mean = FALSE, na.rm = TRUE) {
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]

  w <- if (is.null(w)) rep(1, times = nrow(x)) else as.numeric(w)

  R <- vsum(w * x)
  if (mean) {
    N <- sum(w)
    R <- R / N
  }
  R
}

v_mean <- function(x, w = NULL, na.rm = TRUE) {
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]

  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }

  N <- sum(w)
  xbar <- vsum(x * w) / N
  Rbar <- vlength(xbar)
  xbar / Rbar
}

v_var <- function(x, w = NULL, na.rm = TRUE) {
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]

  Rbar <- vresultant(vnorm(x), w, mean = TRUE, na.rm = FALSE) |>
    vlength()
  1 - Rbar
}

v_sd <- function(x, w = NULL, na.rm = TRUE) {
  Rbar <- vresultant(x, w, mean = TRUE, na.rm = na.rm) |>
    vlength()

  d <- sqrt(log(1 / Rbar^2))
  return(d)
}

v_delta <- function(x, w = NULL, na.rm = TRUE) {
  Rbar <- vresultant(x, w, mean = TRUE, na.rm = na.rm) |> vlength()
  acos(Rbar)
}


v_rdegree <- function(x, w = NULL, na.rm = FALSE) {
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]
  w <- if (is.null(w)) rep(1, times = nrow(x)) else as.numeric(w)

  N <- sum(w)
  Rbar <- vresultant(vnorm(x), w, mean = FALSE, na.rm = FALSE) |> vlength()

  (2 * Rbar - N) / N
}

v_sde <- function(x, w = NULL, na.rm = FALSE) {
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]

  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }

  N <- sum(w)

  if (N < 25) warning("The standard error might not be a good estimator for N < 25")
  xbar <- vsum(x * w) / N
  Rbar <- vlength(xbar)
  mu <- xbar / Rbar

  muv <- mu %*% t(x)

  d <- 1 - sum(muv^2) / N

  angle <- sqrt(d / (N * Rbar^2))
  return(angle)
}

v_confidence_angle <- function(x, w = NULL, alpha = 0.05, na.rm = FALSE) {
  e_alpha <- -log(alpha)
  q <- asin(sqrt(e_alpha) * v_sde(x, w, na.rm = na.rm))
  return(q)
}


#' Statistical estimators of the distribution of a set of vectors
#'
#' @param x object of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param w numeric. Optional weights for each observation.
#' @param alpha numeric. Significance level for the confidence angle (default is 0.05 for a 95% confidence angle).
#' @param na.rm logical. Whether `NA` values should be removed before the computation proceeds.
#' @name stats
#' @details
#' `mean` returns the spherical mean of a set of vectors
#' (object of class of `x`).
#'
#' `var` returns the spherical variance (numeric), based on resultant length
#' (Mardia 1972).
#'
#' `delta` returns the cone angle containing ~63% of the data (in degree if `x` is a `"Plane"` or `"Line"`, or in radians
#' if otherwise). For enough large sample it approaches the angular standard
#' deviation (`"csd"`) of the Fisher statistics.
#'
#' `rdegree` returns the degree of preferred orientation of vectors, range: (0, 1).
#'
#' `sd_error` returns the spherical standard error (numeric). If the number of
#' data is less than 25, if will print a additional message, that the output
#' value might not be a good estimator.
#'
#' `confidence_angle` returns the semi-vertical angle \eqn{q} about the
#' mean \eqn{\mu} (in degree if `x` is a `"Plane"` or `"Line"`, or in radians
#' if otherwise). The \eqn{100(1-\alpha)\%} confidence interval is than given by \eqn{\mu \pm q}.
#' `estimate_k` returns the estimated concentration of the von Mises-Fisher distribution \eqn{\kappa} (after Sra, 2011).#'
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' mean(x)
#' var(x)
#' delta(x)
#' rdegree(x)
#' sd_error(x)
#' confidence_angle(x)
#' estimate_k(x)
#'
#' #' weights:
#' x2 <- Line(c(0, 0), c(0, 90))
#' mean(x2)
#' mean(x2, w = c(1, 2))
#' var(x2)
#' var(x2, w = c(1, 2))
NULL

#' @rdname stats
#' @export
mean.spherical <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))

  x_mean <- Vec3(x) |>
    unclass() |>
    v_mean(w, na.rm) |>
    vec2mat() |>
    Vec3()


  if (is.Line(x)) {
    x_mean <- Line(x_mean)
  } else if (is.Plane(x)) {
    x_mean <- Plane(x_mean)
  }
  rownames(x_mean) <- NULL
  x_mean
}

sd <- function(x, ...) UseMethod("sd")

#' @export
sd.default <- function(x, ...) stats::sd(x, ...)

#' @rdname stats
#' @export
sd.spherical <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  x_sd <- Vec3(x) |>
    unclass() |>
    v_sd(w, na.rm)
  if (is.Line(x) | is.Plane(x)) {
    x_sd * 180 / pi
  } else {
    x_sd
  }
}

var <- function(x, ...) UseMethod("var")

#' @export
var.default <- function(x, ...) stats::var(x, ...)

#' @export
var.data.frame <- function(x, ...) {
  sapply(x, var, ...)
}

#' @rdname stats
#' @export
var.spherical <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  Vec3(x) |>
    unclass() |>
    v_var(w, na.rm)
}

#' @rdname stats
#' @export
confidence_angle <- function(x, w = NULL, alpha = 0.05, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  ca <- Vec3(x) |>
    unclass() |>
    v_confidence_angle(w, alpha, na.rm)

  if (is.Line(x) | is.Plane(x)) {
    ca * 180 / pi
  } else {
    ca
  }
}

#' @rdname stats
#' @export
rdegree <- function(x, w = NULL, na.rm = FALSE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  Vec3(x) |>
    unclass() |>
    v_rdegree(w, na.rm)
}

#' @rdname stats
#' @export
sd_error <- function(x, w = NULL, na.rm = FALSE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  sde <- Vec3(x) |>
    unclass() |>
    v_sde(w, na.rm)

  if (is.Line(x) | is.Plane(x)) {
    sde * 180 / pi
  } else {
    sde
  }
}

#' @rdname stats
#' @export
delta <- function(x, w = NULL, na.rm = TRUE){
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  d <- Vec3(x) |>
    unclass() |>
    v_delta(w, na.rm)

  if (is.Line(x) | is.Plane(x)) {
    d * 180 / pi
  } else {
    d
  }
}

#' @rdname stats
#' @export
estimate_k <- function(x, w = NULL, na.rm = FALSE) {
  Rbar <- vresultant(Vec3(x), w = w, mean = TRUE, na.rm = na.rm) |>
    vlength()

  p <- 3

  (Rbar * (p - Rbar^2)) / (1 - Rbar^2)
}


#' Fisher's statistics
#'
#' Estimates concentration parameter, angular standard deviation, and
#' confidence limit.
#'
#' @inheritParams mean.spherical
#'
#' @returns list, with
#' \describe{
#' \item{`"k"`}{estimated concentration parameter \eqn{\kappa} for the von Mises-Fisher
#' distribution}
#' \item{`"csd"`}{estimated angular standard deviation enclosing 63% of the orientation data. Angle is in degrees if `x` is a spherical object, and raidan if otherwise.}
#' \item{`"a95"`}{Confidence limit for given `p`. Angle is in degrees if `x` is a spherical object, and raidan if otherwise.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' fisher_statistics(x)
fisher_statistics <- function(x, w = NULL, alpha = 0.05, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]
  v <- Vec3(x) |> unclass()

  w <- if (is.null(w)) {
    rep(1, times = nrow(v))
  } else {
    as.numeric(w)
  }

  N <- sum(w)

  R <- vresultant(v, w, na.rm = FALSE) |>
    vlength()

  if (N != R) {
    k <- (N - 1) / (N - R) # fisher's kappa approximation
    csd <- 81 / sqrt(k) # 63%
    csd_95 <- 140 / sqrt(k) # 95%

    term1 <- (N - R) / R
    term2 <- (1 / alpha)^(1 / (N - 1))
    cos_a <- 1 - term1 * (term2 - 1)
    alpha <- acos(cos_a)

    # a95 <- acos(1 - ((N - R) / R) * (20**(1 / (N - 1)) - 1)) # apsg version

    if (is.Line(x) | is.Plane(x)) {
      alpha <- rad2deg(alpha)
      # a95 <- rad2deg(a95)
      # alpha_kn63 <- rad2deg(alpha_kn63)
      # alpha_kn95 <- rad2deg(alpha_kn95)
    } else {
      csd <- deg2rad(csd)
      csd_95 <- deg2rad(csd_95)
    }

    list(k = k, csd = csd, csd_2s = csd_95, alpha = alpha)
  }
}

#' Elliptical concentration and confidence cone estimation
#'
#' @inheritParams mean.spherical
#'
#' @return list
#' \describe{
#'  \item{`k`}{two-column vector containing the estimates for the minimum (\eqn{k_\text{min}}) and maximum concentration (\eqn{k_\text{max}}).}
#'  \item{`a95`}{two-column vector containing the estimates for the minimum and maximum 95% confidence cone.}
#'  \item{`beta`}{The shape factor of the distribution given by the ratio \eqn{\frac{k_\text{min}}{k_\text{max}}}.}
#'  }
#'
#' @export
#'
#' @seealso [inertia_tensor()]
#'
#' @source Borradaile, G. (2003). Spherical-Orientation Data. In: Statistics of
#' Earth Science Data. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-662-05223-5_10
#'
#' @examples
#' set.seed(1234)
#' x <- rfb(100, mu = Line(120, 50), k = 15, A = diag(c(-5, 0, 5)))
#'
#' stereoplot()
#' stereo_point(x)
#'
#' bingham_statistics(x)
bingham_statistics <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  if (isTRUE(na.rm)) x <- x[!rowSums(!is.finite(x)), ]
  v <- Vec3(x)

  w <- if (is.null(w)) {
    rep(1, times = nrow(v))
  } else {
    as.numeric(w)
  }

  n <- sum(w)

  inertia <- inertia_tensor.spherical(v, w)
  abc <- diag(inertia)

  k_ellipse <- numeric(2)
  k_ellipse[2] <- (abc[1] + abc[2]) / (abc[3] + abc[2] - abc[1]) # max
  k_ellipse[1] <- (abc[1] + abc[2]) / (abc[3] - abc[2] + abc[1]) # min

  # theta = x[, 2]
  #
  # k = n / (n - sum(cosd(theta)))
  #
  # A = B = (k * n) / (1+k)
  # C = (2*n) / (1+k)
  #
  # k_ellipse <- c(0, 0)
  # k_ellipse[2] = (A + B) / (C + B - A)
  # k_ellipse[1] = (A + B) / (C - B + A) # min

  a95 <- deg2rad(140) / sqrt(k_ellipse * n)

  if (is.Line(x) | is.Plane(x)) {
    a95 <- rad2deg(a95)
  }

  list(k = k_ellipse, a95 = a95, beta = k_ellipse[1] / k_ellipse[2])
}

#' Test of mean orientations
#'
#' Test against the null-hypothesis that the samples are drawn from the same Fisher population.
#'
#' @param x,y  objects of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @inheritParams mean.spherical
#'
#' @returns list indicating the F-statistic and the p-value.
#'
#' @export
#'
#' @importFrom stats qf
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100, mu = Line(120, 50), k = 20)
#' y <- rvmf(100, mu = Line(180, 45), k = 20)
#'
#' stereoplot()
#' stereo_point(x, col = 1)
#' stereo_point(y, col = 2)
#'
#' fisher_ftest(x, y)
fisher_ftest <- function(x, y, alpha = 0.05, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  if (isTRUE(na.rm)) {
    x <- x[!rowSums(!is.finite(x)), ]
    y <- y[!rowSums(!is.finite(y)), ]
  }

  vx <- Vec3(x)
  vy <- Vec3(y)

  nx <- nrow(vx)
  ny <- nrow(vy)

  Rx <- vresultant(vx) |> vlength()
  Ry <- vresultant(vy) |> vlength()

  R <- vresultant(rbind(vx, vy)) |> vlength()

  stat <- (nx + ny - 2) * ((Rx + Ry - R) / (nx + ny - Rx - Ry))

  df1 <- 2
  df2 <- 2 * (nx + ny - 2)

  crit <- stats::qf(p = alpha, df1 = df1, df2 = df2, lower.tail = FALSE)
  if (stat > crit) {
    message("Reject null-hypothesis")
  } else {
    message("Do not reject null-hypothesis")
  }
  c("F stat" = stat, "p-value" = crit)
}

#' Spherical Linear Interpolation (Slerp)
#'
#' Returns the spherical linear interpolation of points between two vectors
#'
#' @inheritParams fisher_ftest
#' @param t numeric. interpolation factor (`t = [0, 1]`).
#'
#' @note For non-unit vectors the interpolation is not uniform.
#'
#' @details
#' A Slerp path is the spherical geometry equivalent of a path along a line
#' segment in the plane; a great circle is a spherical geodesic.
#'
#' @export
vslerp <- function(x, y, t, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))

  if (isTRUE(na.rm)) {
    x <- x[!rowSums(!is.finite(x)), ]
    y <- y[!rowSums(!is.finite(y)), ]
  }

  vx <- Vec3(x)
  vy <- Vec3(y)

  theta <- angle(x, y)
  slerp <- (x * sin((1 - t) * theta) + y * sin(t * theta)) / sin(theta)

  if (is.Line(x)) {
    Line(slerp)
  } else if (is.Plane(x)) {
    Plane(slerp)
  } else {
    slerp
  }
}
