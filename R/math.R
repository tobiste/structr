#' Vector operations
#'
#' @param x,y,rotaxis numeric vector, array, or object of class `"line"` or `"plane"`
#' @param rotangle numeric. Angle in radians or degrees if `"rotaxis"` is an object of `"line"` or `"plane"`
#' @param w weights.
#' @param l numeric. length to scale `x`
#' @name operations
#' @details
#' 
#' `vlength(x)` returns the length of a vector `x`.
#' 
#' `vsum(x)` returns the vector sum of `x`. If `w` is specified.
#' 
#' `vscale(x)` returns a vector x with length `l`.
#' 
#' `vnorm(x)` returns the normalized vector `x` (i.e. `vlength(x) = 1`)
#' 
#' `vcross(x, y)` returns the cross product of vectors `x` and `y` which is a new vector perpendicular to both.
#' 
#' `vdot(x, y)` returns the dot product (or scalar product) of x and y. `vdot(x, y)` is equivalent to `x %*% t(y)`
#' 
#' `vrotate(x, rotaxis, rotangle)` rotates x about a vector by a specified angle.
#' 
#' `vproject(x,y)` projects x on vector y, i.e. returns a vector that has same orientation as y but length of x.
#' 
#' `vreject(x,y)` Is the vector between `x` and the projected vector of `x` onto `y` (also called the vector resolute of x perpendicular to y)
#' 
#' @examples
#' vec1 <- cbind(1, 0, 0)
#' vec2 <- cbind(0, 0, 1)
#'
#' vlength(vec1)
#' vsum(vec1)
#' vnorm(vec1)
#' vcross(vec1, vec2)
#' vdot(vec1, vec2)
#' vrotate(vec1, vec2, pi / 2)
#' vangle(vec1, vec2)
#' vproject(vec1, vec2)
#' vreject(vec1, vec2)
NULL

#' @rdname operations
#' @export
vlength <- function(x) {
  if (is.spherical(x)) x <- to_vec(x)
  
  sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2) # length of a vector
}

#' @rdname operations
#' @export
vsum <- function(x, w = NULL) {
  if (is.spherical(x)) x <- to_vec(x)
  
  w <- if (is.null(w)) rep(1, times = nrow(x)) else as.numeric(w)
  
  cbind(sum(w * x[, 1]), sum(w * x[, 2]), sum(w * x[, 3])) # vector sum
}

#' @rdname operations
#' @export
vnorm <- function(x) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  }
  xn <- x / vlength(x)
  
  if (transform) to_spherical(xn, class = class) else xn 
}

#' @rdname operations
#' @export
vscale <- function(x, l){
  l/vnorm(x) * x
}

#' @rdname operations
#' @export
vcross <- function(x, y) {
  transform <- is.spherical(x)
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  xxy <- cbind(
    x = x[, 2] * y[, 3] - x[, 3] * y[, 2],
    y = x[, 3] * y[, 1] - x[, 1] * y[, 3],
    z = x[, 1] * y[, 2] - x[, 2] * y[, 1]
  )
  
  if (transform) to_spherical(xxy, class) else xxy
}

#' @rdname operations
#' @export
vdot <- function(x, y) {
  # equivalent to: x %*% t(y)
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) y <- to_vec(y) else {
    y <- vec2mat(y)
  }
  # rowSums(x*y)
  x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
}

#' @rdname operations
#' @export
vrotate <- function(x, rotaxis, rotangle) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(rotaxis)) {
    rotaxis <- to_vec(rotaxis)
    rotangle <- deg2rad(rotangle)
  } else {
    rotaxis <- vec2mat(rotaxis)
  }
  
  rotaxis <- vec2mat(rotaxis) |> vnorm()
  vax <- vcross(rotaxis, x)
  xrot <- x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2 # Helmut
  
  if (transform) {
    to_spherical(xrot, class)
  } else {
    colnames(xrot) <- c("x", "y", "z")
    xrot
  }
}

vrotaxis <- function(x, y) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  xy <- vcross(x, y)
  xy / vnorm(xy)
  
  if (transform) {
    to_spherical(xy, class)
  } else {
    colnames(xy) <- c("x", "y", "z")
    xy
  }
}

# vrotate2 <- function(x, rotaxis, rotangle) {
#   x <- vec2mat(x)
#   rotaxis <- vec2mat(rotaxis)  |> vnorm()
#   vax <- vcross(rotaxis, x)
#   x * cos(rotangle) + vax * sin(rotangle) + rotaxis * vdot(rotaxis, x) *(1-cos(rotangle)) # Rodrigues
# }

#' @rdname operations
#' @export
vangle <- function(x, y) {
  transform <- is.spherical(x)
  if (transform) {
    # classx <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    # classy <- class(y)
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  d <- acos(vdot(x, y))
  
  if (transform) {
    rad2deg(d)
  } else {
    d
  }
}

#' @rdname operations
#' @export
vproject <- function(x, y) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  xpr <- vproject_length(x, y) * y
  
  if (transform) {
    to_spherical(xpr, class)
  } else {
    colnames(xpr) <- c("x", "y", "z")
    xpr
  }
}
vproject_length <- function(x, y) {
  xn <- vnorm(x)
  yn <- vnorm(y)
  # xn * cos(vangle(x, y))
  # as.numeric(xn %*% t(yn))
  vdot(xn, yn)
}

#' @rdname operations
#' @export
vreject <- function(x, y) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  
  x_rej <- x - vproject(x, y)
  
  if (transform) {
    to_spherical(x_rej, class)
  } else {
    colnames(x_rej) <- c("x", "y", "z")
    x_rej
  }
}

#' Orthogonalize two vectors
#'
#' Make two vectors orthogonal by moving them in opposite directions in
#' the plane they span. `vdot(x,y) < 1.0e-4` are considered orthogonal.
#'
#' @param x,y  Coordinates of the two vectors. Any coordinates system
#'
#' @returns  Orthogonalized vectors. In the same coordinate system
#' @export
#'
#' @source https://www.snsn.se/SH/SHcode/Python/BL.py
#'
#' @examples
#' v_orthogonalize(c(1, 0, 0), c(.1, 1, 0))
v_orthogonalize <- function(x, y) {
  # Early exit if already orthogonal
  if (abs(vdot(x, y)) < 1e-4) {
    return(list(x, y))
  }
  
  # Convert to vector matrix form if in spherical
  x_was_sph <- is.spherical(x)
  y_was_sph <- is.spherical(y)
  
  if (x_was_sph) {
    x_class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (y_was_sph) {
    y_class <- class(y)
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  # Compute orthogonalization:
  r <- vnorm(vcross(x, y))  # rotation axis
  a <- vnorm(x)             # first basis
  b <- vcross(r, a)         # second basis
  
  # Project x and y into the rotation coordinate system:
  rv1 <- c(vdot(r, x), vdot(a, x), vdot(b, x))
  rv2 <- c(vdot(r, y), vdot(a, y), vdot(b, y))
  
  # Compute the correction angle:
  ang <- acos(vdot(rv1, rv2))
  ang_correction <- (pi / 2 - ang) / 2
  
  # Rotate in opposite directions around r:
  rw1 <- vrotate(rv1, r, -ang_correction)
  rw2 <- vrotate(rv2, r, ang_correction)
  
  # Transform back to NED:
  w1 <- r * rw1[1] + a * rw1[2] + b * rw1[3]
  w2 <- r * rw2[1] + a * rw2[2] + b * rw2[3]
  
  # Return in the original format
  if (x_was_sph) w1 <- to_spherical(w1, x_class)
  if (y_was_sph) w2 <- to_spherical(w2, y_class)
  
  list(w1, w2)
}


# v_orthogonalize <- function(x, y) {
#   if (abs(vdot(x, y)) < 1.0e-4) {
#     return(list(x, y))
#   } else {
#     transform <- is.spherical(x)
#     if (transform) {
#       classx <- class(x)
#       x <- to_vec(x)
#     } else {
#       x <- vec2mat(x)
#     }
#     if (is.spherical(y)) {
#       classy <- class(y)
#       y <- to_vec(y)
#     } else {
#       y <- vec2mat(y)
#     }
# 
#     # Define the rotation axis (the normal to the v1, v2 plane), use as first
#     # basis vector in the rotation coordinate system.
#     # Remember, the cross product of two non-orthogonal vectors is not a unit
#     # vector
# 
#     r <- vcross(x, y) |> vnorm()
# 
#     a <- vnorm(x)
#     b <- vcross(r, a)
# 
#     # Transform v1 and v2 to the rotation coordinate system.
#     # Transformation matrix is T = (r a b)
#     rv1 <- c(0, 0, 0)
#     rv1[1] <- vdot(r, x)
#     rv1[2] <- vdot(a, x)
#     rv1[3] <- vdot(b, x)
# 
#     rv2 <- c(0, 0, 0)
#     rv2[1] <- vdot(r, y)
#     rv2[2] <- vdot(a, y)
#     rv2[3] <- vdot(b, y)
# 
#     # Rotate rv1 and rv2 half the discrepancy angle in opposite directions
#     # around r
#     ang <- acos(vdot(rv1, rv2))
#     ang <- (pi / 2 - ang) / 2.0
# 
#     rw1 <- vrotate(rv1, r, -ang)
#     rw2 <- vrotate(rv2, r, ang)
# 
#     # Transform the vectors back to NED
#     w1 <- c(0, 0, 0)
#     w1[1] <- r[1] * rw1[1] + a[1] * rw1[2] + b[1] * rw1[3]
#     w1[2] <- r[2] * rw1[1] + a[2] * rw1[2] + b[2] * rw1[3]
#     w1[3] <- r[3] * rw1[1] + a[3] * rw1[2] + b[3] * rw1[3]
# 
#     w2 <- c(0, 0, 0)
#     w2[1] <- r[1] * rw2[1] + a[1] * rw2[2] + b[1] * rw2[3]
#     w2[2] <- r[2] * rw2[1] + a[2] * rw2[2] + b[2] * rw2[3]
#     w2[3] <- r[3] * rw2[1] + a[3] * rw2[2] + b[3] * rw2[3]
# 
#     if (transform) {
#       w1 <- to_spherical(w1, classx)
#       w2 <- to_spherical(w2, classy)
#     }
#     return(list(w1, w2))
#   }
# }

#' Affine transformation of vector by matrix
#'
#' @param x numeric vector, array, or object of class `"line"` or `"plane"`
#' @param A 3x3 matrix
#' @param norm logical. Whether the transformed vector should be normalized
#' (`TRUE`) or not (`FALSE`, the default).
#' @export
#' @examples
#' mat <- cbind(V1 = c(1, 0, 0), V2 = c(0, 1, 0), V3 = c(0, 0, -1))
#' vec <- c(1, 1, 1)
#' vtransform(vec, mat)
vtransform <- function(x, A, norm = FALSE) {
  stopifnot(is.matrix(A))
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  
  xt <- t(A %*% t(x))
  # xt <- t(vdot(A, x))
  if (norm) xt <- vnorm(xt)
  
  if (transform) {
    to_spherical(xt, class)
  } else {
    colnames(xt) <- c("x", "y", "z")
    xt
  }
}


#' Mean resultant of a set of vectors
#'
#' @param x numeric. Can be three element vector or a three column array
#' @param w numerical vector of weights the same length as `x` giving the
#' weights to use for elements of `x`.
#' @param mean logical. Whether the mean resultant (`TRUE`) or resultant
#' (`FALSE`, the default) is returned.
#' @returns if `mean==TRUE`, mean resultant is returned
#' and numeric otherwise.
#' @export
#' @examples
#' x <- rvmf(100, mu = Line(120, 50), k = 5) |> to_vec()
#' vresultant(x, mean = FALSE)
#' vresultant(x, mean = TRUE)
vresultant <- function(x, w = NULL, mean = FALSE) {
  w <- if (is.null(w)) rep(1, times = nrow(x)) else as.numeric(w)
  
  transform <- is.spherical(x)
  if (transform) v <- to_vec(x) else v <- vec2mat(x)
  
  R <- vsum(x, w)
  if (mean) {
    N <- sum(w)
    R <- R / N
  }
  if (transform) to_spherical(R, class(x)) else R
  
}

#' Statistical estimators of the distribution of a set of vectors
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param w numerical vector of weights the same length as `x` giving the
#' weights to use for elements of `x`.
#' @param alpha significance level (0.05 by default)
#' @details
#' `v_mean` returns the spherical mean of a set of vectors
#' (object of class of `x`). The function is shortcut for
#' `vresultant(x, mean = TRUE)` when the argument is a `"plane"` or `"line"` object.
#'
#' `v_var` returns the spherical variance (numeric), based on resultant length
#' (Mardia 1972).
#'
#' `v_delta` returns the cone angle containing ~63% of the data (in degree if `x` is a `"plane"` or `"line"`, or in radians
#' if otherwise). For enough large sample it approaches the angular standard
#' deviation (`"csd"`) of the Fisher statistics.
#'
#' `v_rdegree` returns the degree of preferred orientation of vectors, range: (0, 1).
#'
#' `v_sde` returns the spherical standard error (numeric). If the number of
#' data is less than 25, if will print a additional message, that the output
#' value might not be a good estimator.
#'
#' `v_confidence_angle` returns the semi-vertical angle \eqn{q} about the
#' mean \eqn{\mu} (in degree if `x` is a `"plane"` or `"line"`, or in radians
#' if otherwise). The \eqn{100(1-\alpha)\%} confidence interval is than given by \eqn{\mu \pm q}.
#' `estimate_k` returns the estimated concentration of the von Mises-Fisher distribution \eqn{\kappa} (after Sra, 2011).
#' @seealso [vresultant()], [fisher_statistics()]
#' @name stats
#' @examples
#' set.seed(1234)
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' v_mean(x)
#' v_var(x)
#' v_delta(x)
#' v_rdegree(x)
#' v_sde(x)
#' v_confidence_angle(x)
#' estimate_k(x)
#' fisher_statistics(x)
#'
#' #' weights:
#' x2 <- Line(c(0, 0), c(0, 90))
#' v_mean(x2)
#' v_mean(x2, w = c(1, 2))
#' v_var(x2)
#' v_var(x2, w = c(1, 2))
NULL

#' @rdname stats
#' @export
v_mean <- function(x, w = NULL) {
  transform <- is.spherical(x)
  if (transform) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  
  w <- if (is.null(w)) {
    rep(1, times = nrow(v))
  } else {
    as.numeric(w)
  }
  
  N <- sum(w)
  xbar <- vsum(v, w) / N
  Rbar <- vlength(xbar)
  mu <- xbar / Rbar
  
  if (transform) {
    to_spherical(mu, class(x))
  } else {
    mu
  }
}

#' @rdname stats
#' @export
v_var <- function(x, w = NULL) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(vnorm(v), w, mean = TRUE) |>
    vlength()
  1 - Rbar
}


v_sd <- function(x, w = NULL) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(v, w, mean = TRUE) |>
    vlength()
  
  d <- sqrt(log(1 / Rbar^2))
  if (is.spherical(x)) {
    rad2deg(d)
  } else {
    d
  }
}

#' @rdname stats
#' @export
v_delta <- function(x, w = NULL) {
  transform <- is.spherical(x)
  if (transform) v <- to_vec(x) else v <- vec2mat(x)
  
  Rbar <- vresultant(v, w, mean = TRUE) |> vlength()
  d <- acos(Rbar)
  
  if (transform) rad2deg(d) else d
}

#' @rdname stats
#' @export
v_rdegree <- function(x, w = NULL) {
  v <- if (is.spherical(x)) to_vec(x) else vec2mat(x)
  w <- if (is.null(w)) rep(1, times = nrow(v)) else as.numeric(w)
  
  N <- sum(w)
  Rbar <- vresultant(vnorm(v), w, mean = FALSE) |> vlength()
  
  (2 * Rbar - N) / N
}

#' @rdname stats
#' @export
v_sde <- function(x, w = NULL) {
  v <- if (is.spherical(x)) to_vec(x) else vec2mat(x)
  
  w <- if (is.null(w)) rep(1, times = nrow(v)) else {
    as.numeric(w)
  }
  
  N <- sum(w)
  
  if (N < 25) warning("The standard error might not be a good estimator for N < 25")
  xbar <- vsum(v, w) / N
  Rbar <- vlength(xbar)
  mu <- xbar / Rbar
  
  muv <- mu %*% t(v)
  
  d <- 1 - sum(muv^2) / N
  
  angle <- sqrt(d / (N * Rbar^2))
  if (is.spherical(x)) {
    rad2deg(angle)
  } else {
    angle
  }
}

#' @rdname stats
#' @export
v_confidence_angle <- function(x, w = NULL, alpha = 0.05) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  
  e_alpha <- -log(alpha)
  q <- asin(sqrt(e_alpha) * v_sde(v, w))
  
  if (is.spherical(x)) {
    rad2deg(q)
  } else {
    q
  }
}

v_antipode <- function(x) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  
  # xa <- vrotate(v, c(1, 0, 0), pi)
  xa <- -v
  
  if (transform) {
    to_spherical(xa, class)
  } else {
    colnames(xa) <- c("x", "y", "z")
    xa
  }
}


#' @rdname stats
#' @export
estimate_k <- function(x, w = NULL) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(v, w = w, mean = TRUE) |>
    vlength()
  
  p <- 3
  
  (Rbar * (p - Rbar^2)) / (1 - Rbar^2)
}

#' Fisher's statistics
#'
#' Estimates concentration parameter, angular standard deviation, and
#' confidence limit.
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param w numeric. Weights
#' @param p numeric. Significance level (`0.05` by default), corresponding to
#' \eqn{100 * (1-p) - 95\%} confidence level.
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
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' fisher_statistics(x)
fisher_statistics <- function(x, w = NULL, p = 0.05) {
  transform <- is.spherical(x)
  if (transform) {
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  
  N <- sum(w)
  
  R <- vresultant(x, w) |>
    vlength()
  
  if (N != R) {
    k <- (N - 1) / (N - R) # fisher's kappa approximation
    csd <- 81 / sqrt(k) # 63%
    csd_95 <- 140 / sqrt(k) # 95%
    
    term1 <- (N - R) / R
    term2 <- (1 / p)^(1 / (N - 1))
    cos_a <- 1 - term1 * (term2 - 1)
    alpha <- acos(cos_a)
    
    # a95 <- acos(1 - ((N - R) / R) * (20**(1 / (N - 1)) - 1)) # apsg version
    
    if (transform) {
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
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param w numeric. Weights
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
bingham_statistics <- function(x, w = NULL) {
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  
  n <- sum(w)
  
  inertia <- inertia_tensor(x, w)
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
  
  if (is.spherical(x)) {
    a95 <- rad2deg(a95)
  }
  
  list(k = k_ellipse, a95 = a95, beta = k_ellipse[1] / k_ellipse[2])
}

#' Test of mean orientations
#'
#' Test against the null-hypothesis that the samples are drawn from the same Fisher population.
#'
#' @param x,y  Can be three element vectors, three column arrays, or objects of class `"line"` or `"plane"`
#' @param alpha Significance level
#'
#' @returns list indicating the F-statistic and the p-value.
#'
#' @export
#'
#' @importFrom stats qf
#'
#' @examples
#' set.seed(1234)
#' x <- rvmf(100, mu = Line(120, 50), k = 20)
#' y <- rvmf(100, mu = Line(180, 45), k = 20)
#'
#' ggstereo() +
#'   ggplot2::geom_point(data = gg(x), ggplot2::aes(x, y, color = "x")) +
#'   ggplot2::geom_point(data = gg(y), ggplot2::aes(x, y, color = "y"))
#'
#' fisher_ftest(x, y)
fisher_ftest <- function(x, y, alpha = 0.05) {
  if (is.spherical(x)) {
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  
  nx <- nrow(x)
  ny <- nrow(y)
  
  Rx <- vresultant(x) |> vlength()
  Ry <- vresultant(y) |> vlength()
  
  R <- vresultant(rbind(x, y)) |> vlength()
  
  
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
#' @param x,y numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param t numeric. interpolation factor (`t = [0, 1]`).
#'
#' @note For non-unit vectors the interpolation is not uniform.
#'
#' @details
#' A Slerp path is the spherical geometry equivalent of a path along a line
#' segment in the plane; a great circle is a spherical geodesic.
#'
#' @export
vslerp <- function(x, y, t) {
  transform <- is.spherical(x)
  if (transform) {
    class <- class(x)
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  theta <- vangle(x, y)
  slerp <- (x * sin((1 - t) * theta) + y * sin(t * theta)) / sin(theta)
  
  if (transform) {
    to_spherical(slerp, class)
  } else {
    slerp
  }
}




#' Orientation tensor
#'
#' 3D orientation tensor, which characterize data distribution using
#' eigenvalue method. See (Watson 1966, Scheidegger 1965).
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param norm logical. Whether the tensor should be normalized or not.
#' @param w numeric. weightings
#'
#' @returns matrix
#'
#' @details The normalized orientation tensor is given as \deqn{D = \frac{1}{n} (x_i, y_i, z_i) (x_i, y_i, z_i)^T}
#' n = 1
#'
#' @export
#'
#' @seealso [or_eigen()], [inertia_tensor()]
#'
#' @examples
#' set.seed(1)
#' x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
#' ortensor(x)
ortensor <- function(x, norm = TRUE, w = NULL) {
  if (is.spherical(x)) x <- to_vec(x)
  
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  
  if (norm) {
    # n <- nrow(x)
    n <- sum(w)
  } else {
    n <- 1
  }
  
  or <- (1 / n) * (t(x) %*% x)
  # or <- matrix(nrow = 3, ncol = 3)
  # or[1, 1] <- sum(x[, 1]^2)
  # or[1, 2] <- sum(x[, 1] * x[, 2])
  # or[1, 3] <- sum(x[, 1] * x[, 3])
  # or[2, 1] <- sum(x[, 2] * x[, 1])
  # or[2, 2] <- sum(x[, 2]^2)
  # or[2, 3] <- sum(x[, 2] * x[, 3])
  # or[3, 1] <- sum(x[, 3] * x[, 1])
  # or[3, 2] <- sum(x[, 3] * x[, 2])
  # or[3, 3] <- sum(x[, 3]^2)
  rownames(or) <- colnames(or) <- NULL
  or
}

#' Inertia tensor
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param w numeric. weightings
#'
#' @return 3 x 3 matrix
#' @details \deqn{D = n - (x_i, y_i, z_i) (x_i, y_i, z_i)^T}
#' @export
#' @seealso [ortensor()]
#'
#' @examples
#' set.seed(1)
#' x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
#' inertia_tensor(x)
inertia_tensor <- function(x, w = NULL) {
  if (is.spherical(x)) x <- to_vec(x)
  
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  
  n <- sum(w)
  inertia <- n - (t(x) %*% x)
  
  rownames(inertia) <- colnames(inertia) <- NULL
  inertia
}

#' Eigenvalues and Eigenvectors of a Set of Vectors
#'
#' Decomposition of Orientation Tensor Eigenvectors and Eigenvalues
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param scaled logical. Whether the Eigenvectors should be scaled by the
#' Eigenvalues (only effective if `x` is in Cartesian coordinates).
#'
#' @returns list containing
#' \describe{
#' \item{`values`}{Eigenvalues}
#' \item{`vectors`}{Eigenvectors in coordinate system of `x`}
#' }
#'
#' @seealso [ortensor()], [eigen()]
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' mu <- rvmf(n = 1) |> vec2line()
#' x <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
#' x_eigen <- or_eigen(x)
#' x_eigen
#' stereoplot()
#' stereo_point(x, col = "grey")
#' stereo_point(mu, lab = "mu", col = 4)
#' stereo_point(x_eigen$vectors, col = c(1, 2, 3), lab = c("E1", "E2", "E3"))
or_eigen <- function(x, scaled = FALSE) {
  x_or <- ortensor(x, norm = FALSE)
  x_eigen <- eigen(x_or, symmetric = TRUE)
  x_eigen$vectors <- t(x_eigen$vectors)
  
  if (scaled) {
    x_eigen$vectors[, 1] * x_eigen$values[1]
    x_eigen$vectors[, 2] * x_eigen$values[2]
    x_eigen$vectors[, 3] * x_eigen$values[3]
  }
  
  if (is.line(x)) {
    x_eigen$vectors <- vec2line(x_eigen$vectors)
  } else if (is.plane(x)) {
    x_eigen$vectors <- vec2plane(x_eigen$vectors)
  } else {
    colnames(x_eigen$vectors) <- c("x", "y", "z")
  }
  x_eigen
}

#' Principal Stretches, Strain and Shape Parameters based on the Orientation Tensor.
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#'
#' @importFrom dplyr near
#' @name strain_shape
#'
#' @details \describe{
#' \item{`stretch_ratios`}{Sqrt of eigenvalue ratios}
#' \item{`strain_ratios`}{Log of stretch ratios}
#' \item{`Ramsay`}{strain symmetry (Ramsay, 1983)}
#' \item{`Woodcock`}{Woodcock shape}
#' \item{`Flinn`}{Flinn strain intensity}
#' \item{`Vollmer`}{Point, Girdle, Random, Cylindricity (B), and Uniform Distance (D) Indices (Vollmer 1990; 2020). `D` is a measure of the "distance" from uniformity, and is linear from R to P, and R to G. End members are: uniform D = 0, girdle D = 0.5, cluster D = 1. The 99% level for a test against uniformity for a sample size of 300 is D = 0.1.}
#' \item{`Nadai`}{natural octahedral unit strain and shear (Nadai, 1963)}
#' \item{`Lisle_intensity`}{Intensity index (Lisle, 1985)}
#' \item{`Waterson_intensity`}{strain intensity (Watterson, 1968)}
#' \item{`lode`}{Lode parameter (Lode, 1926)}
#' \item{`kind`}{Descriptive type of ellipsoid}
#' \item{`MAD`}{maximum angular deviation (Kirschvink, 1980)}
#' \item{`US`}{Uniformity statistic of Mardia (1972)}
#' }
#'
#' @seealso [ortensor()], [or_eigen()], [fabric_indexes()]
#'
#' @returns list
#'
#' @references
#' Flinn, Derek.(1963): "On the statistical analysis of fabric diagrams." Geological Journal 3.2: 247-253.
#'
#' Kirschvink, J. (1980): The least-squares line and plane and the analysis of palaeomagnetic data. Geophysical Journal International, 62(3), 699-718.
#'
#' Lisle, Richard J.  (1985): "The use of the orientation tensor for the description and statistical testing of fabrics." Journal of Structural Geology 7.1: 115-117.
#'
#' Lode, Walter (1926): "Versuche über den Einfluß der mittleren Hauptspannung auf das Fließen der Metalle Eisen, Kupfer und Nickel“
#'  (*"Experiments on the influence of the mean principal stress on the flow of the metals iron, copper and nickel"*], Zeitschrift für Physik, vol. 36 (November), pp. 913–939, DOI: 10.1007/BF01400222
#'
#' Mardia, Kantilal Varichand. (1975): "Statistics of directional data." Journal of the Royal Statistical Society Series B: Statistical Methodology 37.3: 349-371.
#'
#' Nadai, A., and Hodge, P. G., Jr. (1963): "Theory of Flow and Fracture of Solids, vol. II." ASME. J. Appl. Mech. December 1963; 30(4): 640. https://doi.org/10.1115/1.3636654
#'
#' Ramsay, John G. (1967): "Folding and fracturing of rocks." Mc Graw Hill Book Company 568.
#'
#' Vollmer, Frederick W. (1990): "An application of eigenvalue methods to structural domain analysis." Geological Society of America Bulletin 102.6: 786-791.
#'
#' Vollmer, Frederick W. (2020): "Representing Progressive Fabric Paths on a Triangular Plot Using a Fabric Density Index and Crystal Axes Eigenvector Barycenters." Geological Society of America Abstracts. Vol. 52.
#'
#' Watterson, Juan. (1968): "Homogeneous deformation of the gneisses of Vesterland, south-west Greenland". No. 78. CA Reitzel.
#'
#' Woodcock, N. H.  (1977): "Specification of fabric shapes using an eigenvalue method." Geological Society of America Bulletin 88.9: 1231-1236.
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 20)
#' principal_stretch(x)
#' principal_strain(x)
#' or_shape_params(x)
NULL

#' @rdname strain_shape
#' @export
principal_stretch <- function(x) {
  # x_eigen <- ortensor(x, norm = TRUE) |> eigen()
  x_eigen <- or_eigen(x)
  s <- sqrt(x_eigen$values)
  names(s) <- c("S1", "S2", "S3")
  return(s)
}

#' @rdname strain_shape
#' @export
principal_strain <- function(x) {
  e <- principal_stretch(x) |> log()
  names(e) <- c("e1", "e2", "e3")
  return(e)
}

#' @rdname strain_shape
#' @export
or_shape_params <- function(x) {
  eig <- ortensor(x, norm = TRUE) |> eigen()
  s <- principal_stretch(x)
  e <- principal_strain(x)
  # names(s) <- names(e) <- NULL
  
  Rxy <- s[1] / s[2]
  Ryz <- s[2] / s[3]
  Rxz <- s[1] / s[3]
  stretch_ratios <- c(Rxy = Rxy, Ryz = Ryz, Rxz = Rxz)
  
  e12 <- e[1] - e[2]
  e13 <- e[1] - e[3]
  e23 <- e[2] - e[3]
  strain_ratios <- c(e12 = e12, e13 = e13, e23 = e23)
  
  shape <- K <- e12 / e23 # strain symmetry (Ramsay, 1983) / Woodcock shape
  
  goct <- 2 * sqrt(e12^2 + e23^2 + e13^2) / 3 # natural octahedral unit shear (Nadai, 1963)
  eoct <- sqrt(3) * goct / 2 # natural octahedral unit strain (Nadai, 1963)
  Nadai <- c(goct = goct, eoct = eoct)
  
  lode <- ifelse((e[1] - e[3]) > 0, (2 * e[2] - e[1] - e[3]) / (e[1] - e[3]), 0)
  
  kind <- dplyr::case_when(
    dplyr::near(eoct, 0) ~ "O",
    lode < -0.75 ~ "L",
    lode > 0.75 ~ "S",
    lode < -0.15 ~ "LLS",
    lode > 0.15 ~ "SSL",
    TRUE ~ "LS"
  )
  
  # Vollmer
  N <- nrow(x)
  P <- eig$values[1] - eig$values[2] #  Point index (Vollmer, 1990)
  G <- 2 * (eig$values[2] - eig$values[3]) #  Girdle index (Vollmer, 1990)
  R <- 3 * eig$values[3] # Random index (Vollmer, 1990)
  B <- P + G #  Cylindricity index (Vollmer, 1990)
  C <- log(eig$values[1] / eig$values[3])
  # I <- 7.5 * ((eig$values[1] / N - (1 / 3))^2 + (eig$values[2] / N - (1 / 3))^2 + (eig$values[3] / N - (1 / 3))^2)
  I <- 7.5 * sum((eig$values / N - 1/3)^2)
  
  # us <- (15 * N / 2) * sum((eig$values[1] - 1 / 3)^2, (eig$values[2] - 1 / 3)^2, (eig$values[3] - 1 / 3)^2) # Uniformity statistic of Mardia
  us <- (15 * N / 2) * sum((eig$values - 1/3)^2)  # Mardia uniformity statistic
  D <- sqrt(us / (5 * N)) # D of Vollmer 2020
  
  Vollmer <- c(P = P, G = G, R = R, B = B, C = C, I = I, D = D)
  
  Lisle_intensity <- 7.5 * sum((eig$values - 1 / 3)^2)
  
  aMAD_l <- atand(sqrt((1 - eig$values[1]) / (eig$values[1]))) # approximate angular deviation from the major axis along E1
  aMAD_p <- atand(sqrt((eig$values[3]) / (1 - eig$values[3]))) # approximate deviation from the plane normal to E3
  aMAD <- ifelse(shape > 1, aMAD_l, aMAD_p)
  
  MAD_l <- atand(sqrt((eig$values[2] + eig$values[3]) / (eig$values[1]))) # Return maximum angular deviation (MAD) of linearly distributed vectors (Kirschvink 1980)
  MAD_p <- atand(sqrt(eig$values[3] / eig$values[2] + eig$values[3] / eig$values[1])) # maximum angular deviation (MAD) of planarly distributed vectors (Kirschvink 1980).
  MAD <- ifelse(shape > 1, MAD_l, MAD_p) #  maximum angular deviation (MAD)
  
  k <- (Rxy - 1) / (Ryz - 1) #  strain symmetry
  d <- sqrt((Rxy - 1)^2 + (Ryz - 1)^2) # strain intensity
  Flinn <- c(intensity = d, symmetry = k)
  
  D <- e12^2 + e23^2 # strain intensity
  Ramsay <- c(intensity = D, symmetry = K)
  Woodcock <- c(strength = e13, shape = K)
  Watterson_intensity <- Rxy + Ryz - 1
  
  # JPF
  # hom.dens <- projection(hom.cpo, upper.proj(hom.cpo), stereonet)
  # # hom.kde <- if(bw != "NA"){kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), h = bw/100*2.4, lims = kde.lims)$z} else{kde2d(unlist(hom.dens[[3]][1]), unlist(hom.dens[[3]][2]), lims = kde.lims)$z}
  # hom.kde <- kde2d(unlist(hom.dens[[1]]), unlist(hom.dens[[2]]), h = bw / 100 * 2, lims = kde.lims)$z
  # hom.norm <- norm(hom.kde, type = "2")
  # dens.norm <- norm(kde, type = "2")
  
  list(
    stretch_ratios = stretch_ratios,
    strain_ratios = strain_ratios,
    Vollmer = Vollmer,
    Flinn = Flinn,
    Ramsay = Ramsay,
    Woodcock = Woodcock,
    Watterson_intensity = Watterson_intensity, # strain intensity (Watterson, 1968)
    Lisle_intensity = Lisle_intensity, # Intensity index (Lisle, 1985).
    Nadai = Nadai,
    Lode = lode, # Lode parameter (Lode, 1926),
    kind = kind, # descriptive type of ellipsoid
    MAD_approx = as.numeric(aMAD), # approximate deviation according to shape
    MAD = as.numeric(MAD), #  maximum angular deviation (MAD)
    US = us
  )
}

# projection <- function(cpo, cart.extra, stereonet = "area") {
#   if (stereonet == "area") {
#     R <- 1
#     x <- R * sqrt(2) * sin(radians(45 - (cpo[, 2] / 2))) * sin(radians(cpo[, 1]))
#     y <- R * sqrt(2) * sin(radians(45 - (cpo[, 2] / 2))) * cos(radians(cpo[, 1]))
#     x2 <- R * sqrt(2) * sin(radians(45 - (cart.extra[, 2] / 2))) * sin(radians(cart.extra[, 1]))
#     y2 <- R * sqrt(2) * sin(radians(45 - (cart.extra[, 2] / 2))) * cos(radians(cart.extra[, 1]))
#   } else {
#     R <- 1
#     x <- R * tan(radians(45 - (cpo[, 2] / 2))) * sin(radians(cpo[, 1]))
#     y <- R * tan(radians(45 - (cpo[, 2] / 2))) * cos(radians(cpo[, 1]))
#     x2 <- R * tan(radians(45 - (cart.extra[, 2] / 2))) * sin(radians(cart.extra[, 1]))
#     y2 <- R * tan(radians(45 - (cart.extra[, 2] / 2))) * cos(radians(cart.extra[, 1]))
#   }
#   all <- data.frame(x = c(x, x2), y = c(y, y2))
#   all$circ <- all$x^2 + all$y^2
#   caxes <- all[-which(all$circ > (1.2)^2), ]
#   caxes <- list(x, y, caxes)
# }


#' Centering vectors
#'
#' Rotate vector object to position that eigenvectors are parallel to
#' axes of coordinate system: E3||X (north-south), E2||X(east-west),
#' E1||X(vertical)
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param max_vertical Whether the maximum of the von Mises-Fisher distribution
#' is already vertical or not.
#'
#' @returns Object of class of `x`
#'
#' @export
#'
#' @seealso [or_eigen()]
#'
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 20)
#' x_centered <- center(x)
#' stereoplot()
#' stereo_point(x, col = "grey")
#' # stereo_point(mu, col = "red")
#' stereo_point(x_centered, col = "black")
#' stereo_point(Line(c(0, 90, 180), c(0, 0, 90)), col = 2:4, lab = c("E3", "E2", "E1"))
center <- function(x, max_vertical = FALSE) {
  transform <- is.spherical(x)
  x_cart <- if (transform) to_vec(x) else vec2mat(x)
  x_eigen <- or_eigen(x_cart)
  
  # x_trans <- matrix(nrow = nrow(x), ncol = 3)
  # for (i in 1:nrow(x)) {
  #   x_trans[i, ] <- vtransform(x_cart[i, ], x_eigen$vectors, norm = TRUE)
  # }
  x_trans <- t(apply(x_cart, 1, vtransform, A = x_eigen$vectors, norm = TRUE))
  
  if (!max_vertical) x_trans <- vrotate(x_trans, cbind(0, -1, 0), pi / 2) 
  if (transform) x_trans <- to_spherical(x_trans, class(x))
  x_trans
}
