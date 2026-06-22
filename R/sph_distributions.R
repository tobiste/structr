#' von Mises-Fisher Distribution
#'
#' Density and random generation for the spherical normal distribution with mean
#' and concentration parameter (\eqn{\kappa}) .
#'
#' @param n integer. number of random samples to be generated
#' @inheritParams sph_mean
#' @param mu Mean vector. object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
#'  rows are the observations and the columns are the coordinates.
#' @param k numeric. The concentration parameter (\eqn{\kappa}) of the von
#' Mises-Fisher distribution
#' @param method character. Algorithm to generate random vectors from a Fisher distribution. 
#' Either `"geologyGeometry"` (the default) to pick the `rayFisher()` algorithm  from the *geologyGeometry* code compilation, or
#' `"rotasym"` to pick the [rotasym::r_vMF()] algorithm from the *rotasym* package.
#'
#' @source Adapted fom [rotasym::r_vMF()] and [rotasym::d_vMF()], and 
#' `geologyGeometry` by Davis, J.R.
#' @importFrom rotasym r_vMF d_vMF
#
#' @name vonmises-fisher
#' @family random
#'
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100, mu = Ray(120, 50), k = 5)
#' contour(x)
#' points(x)
#'
#' dx <- dvmf(x, mu = Ray(120, 50))
#' head(dx)
#'
#' plot(x, col = assign_col(dx))
NULL

#' @rdname vonmises-fisher
#' @export
rvmf <- function(n = 100, mu = Vec3(1, 0, 0), k = 5, method = c('geologyGeometry', 'rotasym')) {
  stopifnot(nrow(mu) == 1)
  method <- match.arg(method)
  
  if(method == 'rotasym'){
  muv <- Vec3(mu) |>
    unclass() |>
    c()
  res <- rotasym::r_vMF(n, muv, k)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)
  } else {
  r <- rayFisher(n = n, mu = as.vector(Vec3(mu)), kappa = k)
  res <- do.call(rbind, r) |> 
    Vec3()
  }

  if (is.Line(mu) | is.Ray(mu)) {
    Ray(res)
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
}


#' Uniformly Distributed Vectors on the Sphere
#'
#' Create uniformly distributed vectors on the sphere
#'
#' @param class character. Coordinate class of the output vectors.
#' @param n integer. number of random samples to be generated
#' @param method character. The algorithm for generating uniformly distributed
#' vectors. Either `"geologyGeometry"` (the default) for generating random points in Cartesian coordinates 
#' (as in the `geologyGeometry` code compilation), 
#' `"sfs"` for the "Spherical Fibonacci Spiral points on a sphere",
#' `"gss"` for "Golden Section Spiral points on a sphere", or the algorithm
#' [rotasym::r_unif_sphere()] from the rotasym package.
#' @returns object of class specified by `"class"` argument
#' @details
#'  `"sfs"` algorithm is from on John Burkardt (http://people.sc.fsu.edu/~jburkardt/),
#'  `"gss` is from  http://www.softimageblog.com/archives/115
#'
#' @importFrom rotasym r_unif_sphere
#' 
#' @source Adapted fom [rotasym::r_unif_sphere()] and `rayUniform()` from 
#' `geologyGeometry` by Davis, J.R.
#'
#' @family random
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x1 <- runif.spherical(n = 100, "Ray", method = "sfs")
#' contour(x1)
#'
#' x2 <- runif.spherical(n = 100, "Ray", method = "gss")
#' contour(x2)
#'
#' x3 <- runif.spherical(n = 100, "Ray", method = "rotasym")
#' contour(x3)
#' 
#' x4 <- runif.spherical(n = 100, "Ray", method = "geologyGeometry")
#' contour(x4)
runif.spherical <- function(n = 100, class = c("Vec3", "Ray", "Line", "Plane"), method = c("geologyGeometry", "gss", "sfs", "rotasym")) {
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
  } else if(method == "rotasym") {
    dc <- rotasym::r_unif_sphere(n, p = 3)
  } else {
    dc <- do.call(rbind, rayUniform(n))
  }

  dc <- Vec3(dc)
  if (class == "Line") {
    Line(dc)
  } else if (class == "Plane") {
    Plane(dc)
  } else if(class == "Ray") {
    Ray(dc)
  } else {
    dc
  }
}

#' Spherical Fisher-Bingham distribution
#'
#' Simulation of random values from a spherical Fisher-Bingham distribution.
#'
#' @inheritParams rvmf
#' @param A symmetric matrix
#' @source Adapted from [Directional::rfb()]
#' @importFrom Directional rfb
#'
#' @family random
#' @export
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
#' 
#' contour(x)
#' points(x)
rfb <- function(n = 100, mu = Vec3(1, 0, 0), k = 5, A) {
  muv <- Vec3(mu) |>
    unclass() |>
    c()

  res <- Directional::rfb(n = n, k = k, m = muv, A = A)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)

  if (is.Line(mu) | is.Ray(mu)) {
    Line(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}

#' Bingham Distribution
#'
#' Simulation of random values from a Bingham distribution with any given symmetric matrix A.
#'
#' @inheritParams rfb
#' @param eigenvalues numeric, three-element vector. Eigenvalues of the diagonal symmetric matrix of the Bingham distribution.
#' @param .class character.
#'
#' @returns  A spherical object of class `.class` and length `n`
#' @name rbing
#'
#' @importFrom Directional rbingham f.rbing
#' @family random
#'
#' @references
#' Fallaize C. J. and Kypraios T. (2016). Exact bayesian inference for the Bingham distribution. Statistics and Computing, 26(1): 349–360. http://arxiv.org/pdf/1401.2894v1.pdf
#'
#' @examples
#' a <- cov(iris[, 1:3])
#' r <- rbingham(100, a, "Line")
#' contour(r)
#' points(r)
#'
#' re <- rbingham_eig(100, c(100, 1, 0), "Line")
#' contour(re)
#' points(re)
NULL

#' @rdname rbing
#' @export
rbingham <- function(n, A, .class = c("Vec3", "Line", "Ray", "Plane")) {
  .class <- match.arg(.class)
  Directional::rbingham(n, A) |>
    as.Vec3() |>
    Spherical(.class)
}

#' @rdname rbing
#' @export
rbingham_eig <- function(n, eigenvalues, .class = c("Vec3", "Line", "Ray", "Plane")) {
  .class <- match.arg(.class)
  ne <- length(eigenvalues)
  stopifnot(ne <= 3)

  if (ne == 3) {
    eigenvalues <- sort(eigenvalues, decreasing = TRUE)
    eigenvalues <- eigenvalues - min(eigenvalues)
    eigenvalues <- eigenvalues[1:2]
  }

  r <- Directional::f.rbing(n, eigenvalues, fast = TRUE)
  Spherical(as.Vec3(r$X), .class)
}


#' Kent Distribution
#'
#' Simulation of random values from a spherical Kent distribution.
#'
#' @inheritParams rvmf
#' @param b numeric. \eqn{\beta} (ellipticity): \eqn{0 \leq \beta < \kappa}
#'
#' @importFrom Directional rkent
#' @family random
#' @export
#' @examples
#' set.seed(20250411)
#' r <- rkent(100, mu = Ray(120, 50), k = 5, b = 1)
#' 
#' contour(r)
#' points(r)
rkent <- function(n, mu = Vec3(1, 0, 0), k = 5, b) {
  muv <- Vec3(mu) |> c()
  res <- Directional::rkent(n = n, k = k, m = muv, b = b)
  colnames(res) <- c("x", "y", "z")
  res <- Vec3(res)

  if (is.Line(mu) | is.Ray(mu)) {
    Ray(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}


#' Watson Distribution
#' 
#' A naive acceptance-rejection sampling algorithm, based on bounding the 
#' density (with respect to the distance from `mu`) with a constant. For large 
#' `kappa`, this method grows inefficient. For `kappa` == 100, about 13 tries 
#' are needed per success. For `kappa` == -100, about 18 tries are needed.
#'
#' @inheritParams rvmf
#'
#' @returns vector of class `mu` of length `n`
#' @family random
#' @export
#' 
#' @source `geologyGeometry` by Davis, J.R.
#'
#' @examples
#' set.seed(20250411)
#' r <- rwatson(100, mu = Ray(120, 50), k = 10)
#' 
#' contour(r)
#' points(r)
rwatson <- function(n, mu, k){
  muv <- Vec3(mu) |> 
    as.vector()
  
  r <- lineWatson(n = n, mu = muv, kappa = k)
  res <- do.call(rbind, r) |> 
    as.Vec3() 
  
  if (is.Line(mu) | is.Ray(mu)) {
    Line(res)
  } else if (is.Plane(mu)) {
    Plane(res)
  } else {
    res
  }
}


#' Random Rotation Matrices
#'
#' Random sample of matrices in SO(3).
#'
#' @inheritParams rvmf
#'
#' @returns list of rotation matrices
#' @export
#'
#' @family random
#'
#' @examples
#' set.seed(20250411)
#' # Generate 10 random SO(3) rotation matrices
#' r <- rrot(10)
#'
#' # convert SO(3) matrices to "Pair"
#' rp <- lapply(r, rot2pair) |> lapply(unclass)
#' rp <- do.call(rbind, args = rp) |>
#'   as.Pair()
#'   
#' # plot pairs
#' plot(rp)
rrot <- function(n) {
  p <- 3
  a <- c(1, numeric(p - 1))
  A <- array(dim = c(p, p, n))
  Ip <- diag(p)
  rotl <- lapply(seq_len(n), function(i) {
    b <- stats::rnorm(p)
    b <- b / sqrt(sum(b^2))
    ca <- a - b * b[1]
    ca <- ca / sqrt(sum(ca^2))
    B <- tcrossprod(b, ca)
    B <- B - t(B)
    theta <- acos(b[1])
    Ip + sin(theta) * B + (cos(theta) - 1) *
      (tcrossprod(b) + tcrossprod(ca))
  })
  if (n == 1) rotl <- list(A)

  return(lapply(rotl, as.Rotation))
}
