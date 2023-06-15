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
#' @seealso [v_unif()] for alternative algorithms to generate uniform
#' distributed samples on a sphere
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
  if (is.spherical(x)) x <- line2vec(x)
  if (is.spherical(mu)) mu <- line2vec(mu) |> c()

  res <- rotasym::d_vMF(x, mu, k)
}


#' Kent distribution



#' Uniformly distributed vectors
#'
#' Create uniformly distributed vectors using the algorithm *Spherical Fibonacci Spiral points on a sphere* algorithm (John Burkardt) or
#'  *Golden Section Spiral points on a sphere*.
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
