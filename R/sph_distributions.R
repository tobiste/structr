#' von Mises-Fisher Distribution
#'
#' Density and random generation for the spherical normal distribution with mean
#' and concentration parameter (\eqn{\kappa}) .
#'
#' @param n integer. number of random samples to be generated
#' @param x,mu object of class `"Vec3"`, `"Line"` or `"Plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa}) of the von
#' Mises-Fisher distribution
#' @seealso [runif.spherical()] for alternative algorithms to generate uniform
#' distributed samples on a sphere, [rkent()] for Kent distribution,
#' [rfb()] for Fisher-Bingham distribution.
#' @source Adapted fom [rotasym::r_vMF()] and [rotasym::d_vMF()]
#' @importFrom rotasym r_vMF d_vMF
# #' @importFrom Directional rvmf dvmf
#' @name vonmises-fisher
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100, mu = Line(120, 50), k = 5)
#' dx <- dvmf(x, mu = Line(120, 50)); head(dx)
#' 
#' plot(x, col = assign_col(dx))
NULL

#' @rdname vonmises-fisher
#' @export
rvmf <- function(n = 100, mu = Vec3(1, 0, 0), k = 5) {
  stopifnot(nrow(mu) == 1)
  muv <- Vec3(mu) |>
    unclass() |>
    c()
  res <- rotasym::r_vMF(n, muv, k)
  # res <- Directional::rvmf(n, muv, k)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)

  if (is.Line(mu)) {
    Line(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}

#' @rdname vonmises-fisher
#' @export
dvmf <- function(x, mu, k = 5) {
  stopifnot(nrow(mu) == 1)
  xv <- Vec3(x) |> unclass()
  muv <- Vec3(mu) |>
    unclass() |>
    c()
  rotasym::d_vMF(xv, muv, k)
  # Directional::dvmf(xv, muv, k) |> as.numeric()
}



#' Uniformly distributed vectors
#'
#' Create uniformly distributed vectors using the algorithm
#' *Spherical Fibonacci Spiral points on a sphere* algorithm (John Burkardt) or
#' *Golden Section Spiral points on a sphere*.
#'
#' @param class character. Coordinate class of the output vectors.
#' @param n integer. number of random samples to be generated
#' @param method character. The algorithm for generating uniformly distributed
#' vectors. Either `"sfs"` for the "Spherical Fibonacci Spiral points on a sphere",
#' `"gss"` for "Golden Section Spiral points on a sphere", or the algorithm
#' [rotasym::r_unif_sphere()] from the rotasym package.
#' @returns object of class specified by `"class"` argument
#' @details
#'  `"sfs"` algorithm is from on John Burkardt (http://people.sc.fsu.edu/~jburkardt/),
#'  `"gss` is from  http://www.softimageblog.com/archives/115
#' @seealso [rvmf()] to draw samples from the von Mises Fisher distribution
#' around a specified mean vector. [rkent()] to draw from a Kent-distribution. [rfb()] to draw from a Fisher-Bingham distribution.
#' @importFrom rotasym r_unif_sphere
#' @export
#' @examples
#' set.seed(20250411)
#' x1 <- runif.spherical(n = 100, "Line", method = "sfs")
#' plot(x1)
#'
#' x2 <- runif.spherical(n = 100, "Line", method = "gss")
#' plot(x2)

#' x3 <- runif.spherical(n = 100, "Line", method = "rotasym")
#' plot(x3)
runif.spherical <- function(n = 100, class = c("Vec3", "Line", "Plane"), method = c("gss", "sfs", "rotasym")) {
  method <- match.arg(method)
  class <- match.arg(class)

  if (method == "sfs") {
    # Spherical Fibonacci Spiral points on a sphere
    phi <- (1 + sqrt(5)) / 2
    i2 <- 2 * (seq(n) - 1) - n + 1
    theta <- 2 * pi * i2 / phi
    sp <- i2 / n
    cp <- sqrt((n + i2) * (n - i2)) / n
    dc <- cbind(x = cp * sin(theta), y = cp * cos(theta), z = sp)
  } else if (method == "sfs") {
    # Golden Section Spiral points on a sphere
    inc <- pi * (3 - sqrt(5))
    off <- 2 / n
    k <- seq(n) - 1
    y <- k * off - 1 + (off / 2)
    r <- sqrt(1 - y * y)
    phi <- k * inc
    dc <- cbind(x = cos(phi) * r, y = y, z = sin(phi) * r)
  } else {
    dc <- rotasym::r_unif_sphere(n, p = 3)
  }

  dc <- Vec3(dc)
  if (class == "Line") {
    Line(dc)
  } else if (class == "Plane") {
    Plane(dc)
  } else {
    dc
  }
}

#' Spherical Fisher-Bingham distribution
#'
#' Simulation of random values from a spherical Fisher-Bingham distribution.
#'
#' @param n integer. number of random samples to be generated
#' @param mu object of class `"Vec3"`, `"Line"` or `"Plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa})
#' @param A symmetric matrix
#' @source Adapted from [Directional::rfb()]
#' @importFrom Directional rfb
#'
#' @seealso [rvmf()] to draw samples from the von Mises Fisher distribution
#' around a specified mean vector. [rkent()] to draw from a Kent-distribution.
#' [runif.spherical()] to draw from a a spherical uniform distribution.
#'
#' @export
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
#' plot(x)
rfb <- function(n = 100, mu = Vec3(1, 0, 0), k = 5, A) {
  muv <- Vec3(mu) |>
    unclass() |>
    c()

  res <- Directional::rfb(n = n, k = k, m = muv, A = A)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)

  if (is.Line(mu)) {
    Line(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}


#' Kent distribution
#'
#' Simulation of random values from a spherical Kent distribution.
#'
#' @inheritParams rvmf
#' @param mu Mean orientation. object of class description `"Vec3()"`, `"Line()"`, or `"Plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa})
#' @param b numeric. \eqn{\beta} (ellipticity): \eqn{0 \leq \beta < \kappa}
#'
#' @source Adapted from [Directional::rkent()]
#' @seealso [rvmf()] to draw samples from the von Mises Fisher distribution
#' around a specified mean vector. [runif.spherical()] to draw from a a spherical uniform distribution.
#'  [rfb()] to draw from a Fisher-Bingham distribution.
#'
#' @importFrom Directional rkent
#' @export
#' @examples
#' set.seed(20250411)
#' x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
#' plot(x)
rkent <- function(n = 100, mu = Vec3(1, 0, 0), k = 5, b) {
  muv <- Vec3(mu) |> c()
  res <- Directional::rkent(n = n, k = k, m = muv, b = b)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)

  if (is.Line(mu)) {
    Line(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}



#' MLE of spherical rotational symmetric distributions
#'
#' Estimates the parameters of a von Mises-Fisher or Kent distribution.
#'
#' @inheritParams ortensor
#' @source Adapted from [Directional::kent.mle()] and [Directional::vmf.mle()]
#' @name dist.mle
#' @importFrom Directional kent.mle vmf.mle
#' @examples
#' x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
#' kent_mle(x)
#' vmf_mle(x)
NULL

#' @rdname dist.mle
#' @export
kent_mle <- function(x) {
  xv <- Vec3(x) |> unclass()
  res <- Directional::kent.mle(xv)
  nm <- colnames(res$G)
  res$G <- t(res$G) |> Vec3()
  rownames(res$G) <- nm
  if (is.Line(x)) {
    res$G <- Line(res$G)
  } else if (is.Plane(x)) {
    res$G <- Plane(res$G)
  }

  res$runtime <- NULL
  return(res)
}

#' @rdname dist.mle
#' @export
vmf_mle <- function(x) {
  xv <- Vec3(x) |> unclass()

  res <- Directional::vmf.mle(xv, fast = TRUE)
  res$mu <- t(res$mu)
  colnames(res$mu) <- c("x", "y", "z")
  res$mu <- Vec3(res$mu)

  if (is.Line(x)) {
    res$mu <- Line(res$mu)
  } else if (is.Plane(x)) {
    res$mu <- Plane(res$mu)
  }

  return(res)
}
