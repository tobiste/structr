#' The von Mises-Fisher Distribution
#'
#' Density and random generation for the spherical normal distribution with mean
#' and kappa.
#'
#' @param n integer. number of random samples to be generated
#' @param x,mu numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa}) of the the von
#' Mises-Fisher distributiuon
#' @name vonmises-fisher
#' @examples
#' # example code
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' stereoplot()
#' stereo_point(x)
NULL

#' @rdname vonmises-fisher
#' @export
rvmf <- function(n = 100, mu = c(1, 0, 0), k = 5) {
  transform <- FALSE
  if (is.structure(mu)) {
    mu <- line2vec(mu) |> c()
    transform <- TRUE
  }

  res <- rotasym::r_vMF(n, mu, k)
  colnames(res) <- c("x", "y", "z")
  if (transform) {
    res |> vec2line()
  } else {
    res
  }
}

#' @rdname vonmises-fisher
#' @export
dvmf <- function(x, mu, k = 5) {
  if (is.structure(x)) x <- line2vec(x)
  if (is.structure(mu)) mu <- line2vec(mu) |> c()

  res <- rotasym::d_vMF(x, mu, k)
}


vec2mat <- function(x) {
  class <- class(x)
  if (is.null(dim(x))) {
    m <- as.matrix(t(x))
  } else {
    m <- as.matrix(x)
  }
  class(m) <- class
  m
}

#' Vector operations
#'
#' @param x,y,rotaxis numeric vector, array, or object of class `"line"` or `"plane"`
#' @param rotangle numeric. Angle in radians or degrees if `"rotaxis"` is an object of `"line"` or `"plane"`
#' @name operations
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
  if (is.structure(x)) {
    x <- to_vec(x)
  }
  sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2) # length of a vector
}

#' @rdname operations
#' @export
vsum <- function(x) {
  if (is.structure(x)) {
    x <- to_vec(x)
  }
  cbind(sum(v[, 1]), sum(v[, 2]), sum(v[, 3])) # vector sum
}

#' @rdname operations
#' @export
vnorm <- function(x) {
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  }
  xn <- x / vlength(x)
  if (transform) {
    to_struct(xn)
  } else {
    xn
  }
}

#' @rdname operations
#' @export
vcross <- function(x, y) {
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }

  vx <- x[, 2] * y[, 3] - x[, 3] * y[, 2]
  vy <- x[, 3] * y[, 1] - x[, 1] * y[, 3]
  vz <- x[, 1] * y[, 2] - x[, 2] * y[, 1]
  xy <- cbind(x = vx, y = vy, z = vz)

  if (transform) {
    to_struct(xy, class)
  } else {
    xy
  }
}

#' @rdname operations
#' @export
vdot <- function(x, y) {
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
}

#' @rdname operations
#' @export
vrotate <- function(x, rotaxis, rotangle) {
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(rotaxis)) {
    rotaxis <- to_vec(rotaxis)
    rotangle <- deg2rad(rotangle)
  } else {
    rotaxis <- vec2mat(rotaxis)
  }

  rotaxis <- vec2mat(rotaxis) |> vnorm()
  vax <- vcross(rotaxis, x)
  xrot <- x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2 # Helmut

  if (transform) {
    to_struct(xrot, class)
  } else {
    colnames(xrot) <- c("x", "y", "z")
    xrot
  }
}

vrotaxis <- function(x, y) {
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }

  xy <- vcross(x, y)
  xy / vnorm(xy)

  if (transform) {
    to_struct(xy, class)
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
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(y)) {
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
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }
  if (is.structure(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }

  xpr <- vproject_length(x, y) * y

  if (transform) {
    to_struct(xpr, class)
  } else {
    colnames(xpr) <- c("x", "y", "z")
    xpr
  }
}
vproject_length <- function(x, y) {
  xn <- vnorm(x)
  # xn * cos(vangle(x, y))
  xn %*% vnorm(y)
}

#' @rdname operations
#' @export
vreject <- function(x, y) {
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }

  x_rej <- x - vproject(x, y)

  if (transform) {
    to_struct(x_rej, class)
  } else {
    colnames(x_rej) <- c("x", "y", "z")
    x_rej
  }
}

#' Affine transformation of vector by matrix
#'
#' @param x numeric vector or array
#' @param A matrix
#' @param norm logical.
#' @examples
#' mat <- cbind(V1 = c(1, 0, 0), V2 = c(0, 1, 0), V3 = c(0, 0, -1))
#' vec <- c(1, 1, 1)
#' vtransform(vec, mat)
vtransform <- function(x, A, norm = TRUE) {
  stopifnot(is.matrix(A))
  transform <- FALSE
  if (is.structure(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  } else {
    x <- vec2mat(x)
  }


  xt <- t(A %*% x)
  if (norm) xt <- vnorm(xt)

  if (transform) {
    to_struct(xt, class)
  } else {
    colnames(xt) <- c("x", "y", "z")
    xt
  }
}


#' Mean resultant of a set of vectors
#'
#' @param x numeric. Can be three element vector or a three column array
#' @param mean logical. Whether the mean resultant (`TRUE`) or resultant
#' (`FALSE`, the default) is returned.
#' @returns if `mean==TRUE`, mean resultant is returned
#' and numeric otherwise.
#' @export
#' @examples
#' x <- rvmf(100, mu = Line(120, 50), k = 5) |> to_vec()
#' vresultant(x, mean = FALSE)
#' vresultant(x, mean = TRUE)
vresultant <- function(x, mean = FALSE) {
  R <- vsum(x)
  if (mean) {
    R <- R / nrow(x)
  }
  R
}

#' Statistical estimators of the distribution of a set of vectors
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
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
#' `v_kappa` returns the estimation for the concentration of the von Mises-Fisher distribution \eqn{\kappa} (after Sra, 2011).
#' @seealso [vresultant()], [fisher_statistics()]
#' @name stats
#' @examples
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' v_mean(x)
#' v_var(x)
#' v_delta(x)
#' v_rdegree(x)
#' v_sde(x)
#' v_confidence_angle(x)
#' v_kappa(x)
NULL

#' @rdname stats
#' @export
v_mean <- function(x) {
  transform <- FALSE
  if (is.structure(x)) {
    v <- to_vec(x)
    transform <- TRUE
  } else {
    v <- vec2mat(x)
  }

  N <- nrow(v)
  xbar <- vsum(v) / N
  Rbar <- vlength(xbar)
  mu <- xbar / Rbar

  if (transform) {
    to_struct(mu, class(x))
  } else {
    mu
  }
}

#' @rdname stats
#' @export
v_var <- function(x) {
  if (is.structure(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(vnorm(v), mean = TRUE) |>
    vlength()
  1 - Rbar
}

v_sd <- function(x) {
  if (is.structure(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }
  Rbar <- vresultant(v, mean = TRUE) |>
    vlength()
  
  sqrt(log(1/Rbar^2))
}



#' @rdname stats
#' @export
v_delta <- function(x) {
  transform <- FALSE
  if (is.structure(x)) {
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
  if (is.structure(x)) {
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
  if (is.structure(x)) {
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
  if (is.structure(x)) {
    v <- to_vec(x)
  } else {
    v <- vec2mat(x)
  }

  e_alpha <- -log(alpha)
  q <- asin(sqrt(e_alpha) * v_sde(v))

  if (is.structure(x)) {
    rad2deg(q)
  } else {
    q
  }
}

#' @rdname stats
#' @export
v_kappa <- function(x) {
  if (is.structure(x)) {
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
  if (is.structure(x)) {
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
  if (is.structure(x)) {
    x <- to_vec(x)
  } else {
    x <- vec2mat(x)
    transform <- TRUE
  }
  if (is.structure(y)) {
    y <- to_vec(y)
  } else {
    y <- vec2mat(y)
  }
  theta <- vangle(x, y)
  slerp <- (x * sin((1 - t) * theta) + y * sin(t * theta)) / sin(theta)

  if (transform) {
    to_struct(slerp, class(x))
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
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
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
#' mu <- Line(120, 50)
#' x <- rvmf(100, mu = mu, k = 5)
#' x_eigen <- or_eigen(x)
#' x_eigen
#' stereoplot()
#' stereo_point(x, col = "grey")
#' stereo_point(mu, lab = "mu")
#' stereo_point(x_eigen$vectors, col = c(1, 2, 3), lab = c("E1", "E2", "E3"))
or_eigen <- function(x, scaled = FALSE) {
  x_or <- ortensor(x, norm = FALSE)
  x_eigen <- eigen(x_or)
  x_eigen$vectors <- t(x_eigen$vectors)

  if (is.line(x)) {
    x_eigen$vectors <- x_eigen$vectors |> vec2line()
  } else if (is.plane(x)) {
    x_eigen$plane <- x_eigen$vectors |> vec2plane()
  } else {
    if (scaled) {
      x_eigen$vectors[, 1] * x_eigen$values[1]
      x_eigen$vectors[, 2] * x_eigen$values[2]
      x_eigen$vectors[, 3] * x_eigen$values[3]
    }
    colnames(x_eigen$vectors) <- c("x", "y", "z")
  }
  x_eigen
}

#' Centering vectors
#'
#' Rotate vector object to position that eigenvectors are parallel to
#' axes of coordinate system: E1||X (north-south), E2||X(east-west),
#' E3||X(vertical)
#'
#' @param x numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param max_vertical Whether the maximum of the von Mises-Fisher distribution
#' is already vertical or not.
#' @returns Object of class of `x`
#' @export
#' @seealso [or_eigen()]
#' @examples
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' x_centered <- center(x)
#' stereoplot()
#' stereo_point(x, col = "grey")
#' stereo_point(x_centered, col = "black")
center <- function(x, max_vertical = FALSE) {
  x_cart <- to_vec(x)
  x_or <- ortensor(x_cart)
  # x_svd <- svd(x_or) # Singular Value Decomposition of a Matrix
  x_eigen <- eigen(x_or)
  x_eigen$vectors <- t(x_eigen$vectors)

  x_trans <- matrix(nrow = nrow(x), ncol = 3)
  for (i in 1:nrow(x)) {
    x_trans[i, ] <- vtransform(x_cart[i, ], x_eigen$vectors)
  }

  if (!max_vertical) {
    x_cent <- vrotate(x_trans, cbind(0, -1, 0), pi / 2)
  } else {
    x_cent <- x_trans
  }
  if (is.structure(x)) {
    x_cent <- to_struct(x_cent, class(x))
  }
  x_cent
}
