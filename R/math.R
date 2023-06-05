DEG2RAD <- function() {
  pi / 180
}

rvmf <- function(n, mu, k, cart = TRUE){
  if(cart){
    Directional::rvmf(n, mu = mu, k = k) |> Directional::euclid()
  } else {
    Directional::rvmf(n, mu = line2vec(mu) |> t(), k = k) |> 
      Directional::euclid() |>
      vec2line()
  }
}



vec2mat <- function(x) {
  if (is.null(dim(x))) {
    as.matrix(t(x))
  } else {
    as.matrix(x)
  }
}


#' length of a vector
vlength <- function(x) {
  sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2)
}

# normalized vector
vnorm <- function(x) {
  x / vlength(x)
}

# cross product of two vectors
vcross <- function(x, y) {
  vx <- x[, 2] * y[, 3] - x[, 3] * y[, 2]
  vy <- x[, 3] * y[, 1] - x[, 1] * y[, 3]
  vz <- x[, 1] * y[, 2] - x[, 2] * y[, 1]
  cbind(x = vx, y = vy, z = vz)
}

# dot product
vdot <- function(x, y) {
  x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
  # x %*% t(y)
}

# rotate vector around a given axis about a given angle (in radians)
vrotate <- function(x, rotaxis, rotangle) {
  vax <- vcross(rotaxis, x)
  x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2
}

# angle between two vectors in radians
vangle <- function(x, y) {
  acos(vnorm(x) %*% vnorm(y))
}

# projection of x on y
vproject <- function(x, y) {
  vproject_length(x, y) * y
}
vproject_length <- function(x, y) {
  xn <- vnorm(x)
  # xn * cos(vangle(x, y))
  xn * (xn %*% vnorm(y))
}

#   Return vector rejection on the vector other
vreject <- function(x, y) {
  x - vproject(x, y)
}

#' mean resultant of a set of vectors
#' @examples
#' x <- rvmf(5, mu = line2vec(c(90, 10)) |> c(), k = 2, cart = TRUE)
#' vresultant(x, mean = TRUE)  
vresultant <- function(x, mean = FALSE) {
  R <- sum(x)
  if (mean) R <- R / length(x)
  R
}

#' Fisher's statistics
#
#         fisher_statistics returns dictionary with keys:
#             `k`    estimated precision parameter,
#             `csd`  estimated angular standard deviation
#             `a95`  confidence limit
#'
#' @examples
#' x <- rvmf(5, mu = line2vec(c(90, 10)) |> c(), k = 2, cart = TRUE)
#' v_fisher_statistics(x)        
v_fisher_statistics <- function(x) {
  N <- length(x)
  R <- abs(vresultant(vnorm(x)))
  if (N != R) {
    k <- (N - 1) / (N - R)
    csd <- 81 / sqrt(k)
    a95 <- acos(1 - ((N - R) / R) * (20**(1 / (N - 1)) - 1))
    list(k = k, csd = csd, a95 = a95)
  }
}

# """Spherical variance based on resultant length (Mardia 1972).
#
#         var = 1 - abs(R) / n
#   
#' @examples
#' x <- rvmf(5, mu = line2vec(c(90, 10)) |> c(), k = 2, cart = TRUE)
#' v_var(x)  
v_var <- function(x) {
  1 - abs(vresultant(vnorm(x), mean = TRUE))
}

# """Cone angle containing ~63% of the data in degrees.
#
#         For enough large sample it approach angular standard deviation (csd)
#         of Fisher statistics
#         """
#' x <- rvmf(5, mu = line2vec(c(90, 10)) |> c(), k = 2, cart = TRUE)
#' v_delta(x)  
v_delta <- function(x) {
  acos(abs(vresultant(x, mean = TRUE)))
}

# """Degree of preferred orientation of vectors in ``FeatureSet``.
#
#         D = 100 * (2 * abs(R) - n) / n
#         """
#' @examples
#' x <- rvmf(5, mu = line2vec(c(90, 10)) |> c(), k = 2, cart = TRUE)
#' v_rdegree(x)
v_rdegree <- function(x) {
  N <- length(x)
  100 * (2 * abs(vresultant(vnorm(x)) - N)) / N
}

# Return a spherical linear interpolation between self and other vector
#
# Note that for non-unit vectors the interpolation is not uniform
vslerp <- function(x, y, t) {
  theta <- vangle(x, yb)
  x * (x * sin((1 - t) * theta) + y * sin(t * theta)) / sin(theta)
}
