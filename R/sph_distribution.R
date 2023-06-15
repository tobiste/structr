#' von Mises-Fisher Distribution
#'
#' Density and random generation for the spherical normal distribution with mean
#' and kappa.
#'
#' @param n integer. number of random samples to be generated
#' @param x,mu numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa}) of the the von
#' Mises-Fisher distributiuon
#' @seealso [v_unif()] for alternative algorithms to generate uniform
#' distributed samples on a sphere, [rkent()] for Kent distribution.
#' @source Adapted fom [rotasym::r_vMF()] and [rotasym::d_vMF()]
#' @importFrom rotasym r_vMF d_vMF
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
  if (is.spherical(mu)) {
    class = class(mu)
    mu <- line2vec(mu) |> c()
    transform <- TRUE
  }

  res <- rotasym::r_vMF(n, mu, k)
  colnames(res) <- c("x", "y", "z")
  if (transform) {
    res |> to_spherical(class)
  } else {
    res
  }
}

#' @rdname vonmises-fisher
#' @export
dvmf <- function(x, mu, k = 5) {
  if (is.spherical(x)) x <- line2vec(x)
  if (is.spherical(mu)) mu <- line2vec(mu) |> c()

  res <- rotasym::d_vMF(x, mu, k)
}



#' Uniformly distributed vectors
#'
#' Create uniformly distributed vectors using the algorithm *Spherical Fibonacci Spiral points on a sphere* algorithm (John Burkardt) or
#' *Golden Section Spiral points on a sphere*.
#'
#' @param class character. Coordinate class of the output vectors.
#' @param n number of observations
#' @param method chaarcter. The algorithm for generating uniformly distributed
#' vectors. Either `"sfs"` for the "Spherical Fibonacci Spiral points on a sphere"
#' or `"gss"` for "Golden Section Spiral points on a sphere".
#' @returns object of class specified by `"class"` argument
#' @details
#'  `"sfs"` algorithm is from on John Burkardt (http://people.sc.fsu.edu/~jburkardt/),
#'  `"gss` is from  http://www.softimageblog.com/archives/115
#' @seealso [rvmf()] to draw samples from the von Mises Fisher distribution
#' around a specified mean vector.
#' @export
#' @examples
#' v_unif("line", n = 100, method = "sfs") |>
#'   ortensor() |>
#'   or_eigen()
#' v_unif("line", n = 100, method = "gss") |>
#'   ortensor() |>
#'   or_eigen()
v_unif <- function(class = NULL, n = 100, method = c("gss", "sfs")) {
  method <- match.arg(method)

  if (method == "sfs") {
    # Spherical Fibonacci Spiral points on a sphere
    phi <- (1 + sqrt(5)) / 2
    i2 <- 2 * (seq(n) - 1) - n + 1
    theta <- 2 * pi * i2 / phi
    sp <- i2 / n
    cp <- sqrt((n + i2) * (n - i2)) / n
    dc <- cbind(x = cp * sin(theta), y = cp * cos(theta), z = sp)
  } else {
    # Golden Section Spiral points on a sphere
    inc <- pi * (3 - sqrt(5))
    off <- 2 / n
    k <- seq(n) - 1
    y <- k * off - 1 + (off / 2)
    r <- sqrt(1 - y * y)
    phi <- k * inc
    dc <- cbind(x = cos(phi) * r, y = y, z = sin(phi) * r)
  }
  if (!is.null(class)) {
    dc <- to_spherical(dc, class)
  }
  dc
}

#' Spherical Fisher-Bingham distribution
#' 
#' Simulation of random values from a spherical Fisher-Bingham distribution.
#'
#' @param n integer. number of random samples to be generated
#' @param mu numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa})
#' @param A symmetric matrix
#' @source Adapted from [Directional::rfb()]
#' @importFrom Directional rfb
#' @examples
#' \dontrun{
#' x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
#' stereoplot()
#' stereo_point(x)
#' }
rfb <- function(n = 100, mu = c(1, 0, 0), k = 5, A){
  transform <- FALSE
  if (is.spherical(mu)) {
    class <- class(mu)
    mu <- line2vec(mu) |> c()
    transform <- TRUE
  }
  res <- Directional::rfb(n = n, k = k, m = mu, A=A)
  colnames(res) <- c("x", "y", "z")
  if (transform) {
    res |> to_spherical(class)
  } else {
    res
  }
}


#' Kent distribution
#' 
#' Simulation of random values from a spherical Kent distribution.
#'
#' @param n integer. number of random samples to be generated
#' @param mu numeric. Can be three element vector, three column array, or an
#' object of class `"line"` or `"plane"`
#' @param k numeric. The concentration parameter (\eqn{\kappa})
#' @param b numeric. \eqn{\beta} (ellipticity): \eqn{0 \leq \beta < \kappa}
#' @source Adapted from [Directional::rkent()]
#' @seealso [rvmf()] for von Mises-Fisher distribution
#' @importFrom Directional rkent
#' @export
#' @examples
#' x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
#' stereoplot()
#' stereo_point(x)
rkent <- function(n = 100, mu = c(1, 0, 0), k = 5, b){
  transform <- FALSE
  if (is.spherical(mu)) {
    class <- class(mu)
    mu <- line2vec(mu) |> c()
    transform <- TRUE
  }
  
  res <- Directional::rkent(n = n, k = k, m = mu, b=b)
  colnames(res) <- c("x", "y", "z")
  if (transform) {
    res |> to_spherical(class)
  } else {
    res
  }
}


#' MLE of spherical rotational symmetric distributions
#' 
#' Estimates the parameters of a von Mises-Fisher or Kent distribution. 
#'
#' @param x numeric. Can be three element vector, three column array, or an 
#' object of class `"line"` or `"plane"`
#' @source Adapted from [Directional::kent.mle()] and [Directional::vmf.mle()]
#' @name dist.mle
#' @importFrom Directional kent.mle vmf.mle
#' @examples
#' x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
#' kent.mle(x)
#' vmf.mle(x)
NULL

#' @rdname dist.mle
#' @export
kent.mle <- function(x){
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  }
  res = Directional::kent.mle(x)
  
  res$G <- t(res$G)
  if(transform){
    res$G <- to_spherical(res$G, class)
  } 
  res$runtime <- NULL
  return(res)
}

#' @rdname dist.mle
#' @export
vmf.mle <- function(x){
  transform <- FALSE
  if (is.spherical(x)) {
    class <- class(x)
    x <- to_vec(x)
    transform <- TRUE
  }
  res = Directional::vmf.mle(x, fast = TRUE)
  
  if(transform){
    res$mu <- to_spherical(res$mu, class)
  } 
  return(res)
}