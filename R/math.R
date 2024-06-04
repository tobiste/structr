#' Vector operations
#'
#' @param x,y,rotaxis numeric vector, array, or object of class `"line"` or `"plane"`
#' @param rotangle numeric. Angle in radians or degrees if `"rotaxis"` is an object of `"line"` or `"plane"`
#' @name operations
#' @details
#' `vdot(x, y)` is equivalent to `x %*% t(y)`
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
  if (is.spherical(x)) {
    x <- to_vec(x)
  }
  sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2) # length of a vector
}

#' @rdname operations
#' @export
vsum <- function(x, w = NULL) {
  if (is.spherical(x)) {
    x <- to_vec(x)
  }
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  cbind(sum(w * x[, 1]), sum(w * x[, 2]), sum(w * x[, 3])) # vector sum
}

#' @rdname operations
#' @export
vnorm <- function(x) {
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  }
  xn <- x / vlength(x)
  if (transform) {
    to_spherical(xn)
  } else {
    xn
  }
}

#' @rdname operations
#' @export
vcross <- function(x, y) {
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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

  if (transform) {
    to_spherical(xxy, class)
  } else {
    xxy
  }
}

#' @rdname operations
#' @export
vdot <- function(x, y) {
  # equivalent to: x %*% t(y)
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    # transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  # rowSums(x*y)
  x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
}

#' @rdname operations
#' @export
vrotate <- function(x, rotaxis, rotangle) {
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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
  transform <- FALSE
  if (is.spherical(x)) {
    classx <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.spherical(y)) {
    classy <- class(y)
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
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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
  if (abs(vdot(x, y)) < 1.0e-4) {
    return(list(x, y))
  } else {
    transform <- FALSE
    if (is.spherical(x)) {
      class <- class(x)
      x <- to_vec(x)
      transform <- TRUE
    } else {
      x <- vec2mat(x)
    }
    if (is.spherical(y)) {
      y <- to_vec(y)
    } else {
      y <- vec2mat(y)
    }

    # Define the rotation axis (the normal to the v1, v2 plane), use as first
    # basis vector in the rotation coordinate system.
    # Remember, the cross product of two non-orthogonal vectors is not a unit
    # vector

    r <- vcross(x, y) |> vnorm()

    a <- vnorm(x)
    b <- vcross(r, a)

    # Transform v1 and v2 to the rotation coordinate system.
    # Transformation matrix is T = (r a b)
    rv1 <- c(0, 0, 0)
    rv1[1] <- vdot(r, x)
    rv1[2] <- vdot(a, x)
    rv1[3] <- vdot(b, x)

    rv2 <- c(0, 0, 0)
    rv2[1] <- vdot(r, y)
    rv2[2] <- vdot(a, y)
    rv2[3] <- vdot(b, y)

    # Rotate rv1 and rv2 half the discrepancy angle in opposite directions
    # around r
    ang <- acos(vdot(rv1, rv2))
    ang <- (pi / 2 - ang) / 2.0

    rw1 <- vrotate(rv1, r, -ang)
    rw2 <- vrotate(rv2, r, ang)

    # Transform the vectors back to NED
    w1 <- c(0, 0, 0)
    w1[1] <- r[1] * rw1[1] + a[1] * rw1[2] + b[1] * rw1[3]
    w1[2] <- r[2] * rw1[1] + a[2] * rw1[2] + b[2] * rw1[3]
    w1[3] <- r[3] * rw1[1] + a[3] * rw1[2] + b[3] * rw1[3]

    w2 <- c(0, 0, 0)
    w2[1] <- r[1] * rw2[1] + a[1] * rw2[2] + b[1] * rw2[3]
    w2[2] <- r[2] * rw2[1] + a[2] * rw2[2] + b[2] * rw2[3]
    w2[3] <- r[3] * rw2[1] + a[3] * rw2[2] + b[3] * rw2[3]

    if (transform) {
      w1 <- to_spherical(w1, classx)
      w2 <- to_spherical(w2, classy)
    }
    return(list(w1, w2))
  }
}


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
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
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
  w <- if (is.null(w)) {
    rep(1, times = nrow(x))
  } else {
    as.numeric(w)
  }
  
  R <- vsum(x, w)
  if (mean) {
    N <- sum(w)
    R <- R / N
  }
  R
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
#' `v_rdegree` returns the degree of preferred orientation of vectors (in %).
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
#' x <- rvmf(100, mu = Line(120, 50), k = 5)

#' v_var(x)
#' v_delta(x)
#' v_rdegree(x)
#' v_sde(x)
#' v_confidence_angle(x)
#' estimate_k(x)
#' 
#' #' weights:
#' x2 <-  Line(c(0, 0), c(0, 90))
#' v_mean(x2)
#' v_mean(x2, w = c(1, 2))
#' v_var(x2)
#' v_var(x2, w = c(1, 2))
NULL

#' @rdname stats
#' @export
v_mean <- function(x, w = NULL) {
  transform <- FALSE
  if (is.spherical(x)) {
    v <- to_vec(x)
    transform <- TRUE
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

  sqrt(log(1 / Rbar^2))
}

#' @rdname stats
#' @export
v_delta <- function(x, w = NULL) {
  transform <- FALSE
  if (is.spherical(x)) {
    v <- to_vec(x)
    transform <- TRUE
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(v, mean = TRUE) |>
    vlength()
  d <- acos(Rbar)

  if (transform) {
    rad2deg(d)
  } else {
    d
  }
}

#' @rdname stats
#' @export
v_rdegree <- function(x) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  N <- nrow(v)
  Rbar <- vresultant(vnorm(v), mean = FALSE) |>
    vlength()

  100 * (2 * Rbar - N) / N
}

#' @rdname stats
#' @export
v_sde <- function(x) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  N <- nrow(v)
  if (N < 25) warning("The standard error might not be a good estimator for N < 25")
  xbar <- vsum(v) / N
  Rbar <- vlength(xbar)
  mu <- xbar / Rbar

  muv <- mu %*% t(v)

  d <- 1 - sum(muv^2) / N

  sqrt(d / (N * Rbar^2))
}

#' @rdname stats
#' @export
v_confidence_angle <- function(x, alpha = 0.05) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }

  e_alpha <- -log(alpha)
  q <- asin(sqrt(e_alpha) * v_sde(v))

  if (is.spherical(x)) {
    rad2deg(q)
  } else {
    q
  }
}



#' @rdname stats
#' @export
estimate_k <- function(x) {
  if (is.spherical(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  # N <- nrow(v)
  # xbar <- vsum(v) / N
  # Rbar <- vlength(xbar)
  Rbar <- vresultant(v, mean = TRUE) |>
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
#' @returns list, with
#' \describe{
#' \item{`"k"`}{estimated concentration parameter \eqn{\kappa} for the von Mises-Fisher
#' distribution}
#' \item{`"csd"`}{estimated angular standard deviation}
#' \item{`"a95"`}{confidence limit}
#' }
#' @examples
#' \dontrun{
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' fisher_statistics(x)
#' }
fisher_statistics <- function(x) {
  transform <- FALSE
  if (is.spherical(x)) {
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  N <- nrow(x)
  R <- x |>
    vresultant() |>
    vlength()

  if (N != R) {
    k <- (N - 1) / (N - R)
    csd <- 81 / sqrt(k)
    a95 <- acos(1 - ((N - R) / R) * (20**(1 / (N - 1)) - 1)) * ifelse(transform, 1 / DEG2RAD(), 1)
    list(k = k, csd = csd, a95 = a95)
  }
}




#' Spherical Linear Interpolation (Slerp)
#'
#' Returns the spherical linear interpolation of points between two vectors
#'
#' @param x,y numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param t numeric. interpolation factor (`t = [0, 1]`).
#' @note For non-unit vectors the interpolation is not uniform
#' @details
#' A Slerp path is the spherical geometry equivalent of a path along a line
#' segment in the plane; a great circle is a spherical geodesic. #'
#' @export
vslerp <- function(x, y, t) {
  transform <- FALSE
  if (is.spherical(x)) {
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
    transform <- TRUE
  }
  if (is.spherical(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  theta <- vangle(x, y)
  slerp <- (x * sin((1 - t) * theta) + y * sin(t * theta)) / sin(theta)

  if (transform) {
    to_spherical(slerp, class(x))
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
#' @returns matrix
#' @export
#' @seealso [or_eigen()]
#' @examples
#' set.seed(1)
#' x <- rvmf(100, mu = Line(120, 50), k = 20)
#' ortensor(x)
ortensor <- function(x, norm = TRUE) {
  if (is.line(x) | is.plane(x)) x <- to_vec(x)
  or <- t(x) %*% x
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
  if (norm) {
    or / nrow(x)
  } else {
    or
  }
}


#' Eigenvalues and eigenvectors of a set of vectors
#'
#' Decomposition of Orientation tensor
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param scaled logical. Whether the eigenvectors should be scaled by the
#' eigenvalues (only effective if `x` is in Cartesian coordinates)
#' @returns list containing
#' \describe{
#' \item{`values`}{Eigenvalues}
#' \item{`vectors`}{Eigenvectors in coordinates system of `x`}
#' }
#' @seealso [ortensor()]
#' @export
#' @examples
#' set.seed(1)
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 20)
#' x_eigen <- or_eigen(x)
#' x_eigen
#' stereoplot()
#' stereo_point(x, col = "grey")
#' stereo_point(mu, lab = "mu", col = 4)
#' stereo_point(x_eigen$vectors, col = c(1, 2, 3), lab = c("E1", "E2", "E3"))
or_eigen <- function(x, scaled = FALSE) {
  x_or <- ortensor(x, norm = TRUE)
  # x_eigen <- eigen(x_or)
  # x_eigen$vectors <- t(x_eigen$vectors)

  x_svd <- svd(x_or) # Singular Value Decomposition of a Matrix
  x_eigen <- list(
    values = x_svd$d,
    vectors = t(x_svd$u)
  )


  if (scaled) {
    x_eigen$vectors[, 1] * x_eigen$values[1]
    x_eigen$vectors[, 2] * x_eigen$values[2]
    x_eigen$vectors[, 3] * x_eigen$values[3]
  }

  if (is.line(x)) {
    x_eigen$vectors <- x_eigen$vectors |> vec2line()
  } else if (is.plane(x)) {
    x_eigen$vectors <- x_eigen$vectors |> vec2plane()
  } else {
    colnames(x_eigen$vectors) <- c("x", "y", "z")
  }
  x_eigen
}

#' Principal stretches, strain and shape parameters based on the orientation tensor.
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#'
#' @importFrom dplyr near
#' @name strain_shape
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
  names(s) <- names(e) <- NULL

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

  if (dplyr::near(eoct, 0)) {
    kind <- "O"
  } else if (lode < -0.75) {
    kind <- "L"
  } else if (lode > 0.75) {
    kind <- "S"
  } else if (lode < -0.15) {
    kind <- "LLS"
  } else if (lode > 0.15) {
    kind <- "SSL"
  } else {
    kind <- "LS"
  }

  # Vollmer
  P <- eig$values[1] - eig$values[2] #  Point index (Vollmer, 1990)
  G <- 2 * (eig$values[2] - eig$values[3]) #  Girdle index (Vollmer, 1990)
  R <- 3 * eig$values[3] # Random index (Vollmer, 1990)
  B <- P + G #  Cylindricity index (Vollmer, 1990)

  Vollmer <- c(P = P, G = G, R = R, B = B)

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
    MAD = as.numeric(MAD) #  maximum angular deviation (MAD)
  )
}

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
#' @returns Object of class of `x`
#' @export
#' @seealso [or_eigen()]
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
  x_cart <- to_vec(x)
  x_or <- ortensor(x_cart)
  x_svd <- svd(x_or) # Singular Value Decomposition of a Matrix
  x_eigen <- list(
    values = x_svd$d,
    vectors = t(x_svd$u)
  )
  # x_eigen <- eigen(x_or)
  # x_eigen$vectors <- t(x_eigen$vectors)

  x_trans <- matrix(nrow = nrow(x), ncol = 3)
  for (i in 1:nrow(x)) {
    x_trans[i, ] <- vtransform(x_cart[i, ], x_eigen$vectors, norm = TRUE)
  }

  if (!max_vertical) {
    x_cent <- vrotate(x_trans, cbind(0, -1, 0), pi / 2)
  } else {
    x_cent <- x_trans
  }
  if (is.spherical(x)) {
    x_cent <- to_spherical(x_cent, class(x))
  }
  x_cent
}
