#' @include vector_math.R
NULL

Test <- function(name, statistic, p_value = NULL, df = NULL,
                 critical = NULL, level = 0.95, source = "",
                 extra = list()) {
  # Result of a hypothesis test. p_value is NULL when the test only has
  # critical values; extra holds further statistics (e.g. resultant length R).
  structure(
    list(
      name      = name,
      statistic = statistic,
      p_value   = p_value,
      df        = df,
      critical  = critical,
      level     = level,
      source    = source,
      extra     = extra
    ),
    class = "Test"
  )
}

test_reject <- function(x, level = NULL, ...) {
  stopifnot(inherits(x, 'Test'))
  # TRUE if the null hypothesis is rejected at confidence `level`
  # (default: the level of the test).
  if (is.null(level)) level <- x$level
  if (!is.null(x$p_value)) {
    return(x$p_value < 1 - level)
  }
  if (level != x$level) {
    stop(sprintf("%s has a critical value for level %s only", x$name, x$level),
         call. = FALSE)
  }
  if (is.null(x$critical)) {
    stop(sprintf("%s has neither a p-value nor a critical value", x$name),
         call. = FALSE)
  }
  x$statistic > x$critical
}

#' @exportS3Method base::as.list
as.list.Test <- function(x, ...) {
  c(
    list(
      name      = x$name,
      statistic = x$statistic,
      p_value   = x$p_value,
      df        = x$df,
      critical  = x$critical,
      level     = x$level,
      reject    = test_reject(x)
    ),
    x$extra
  )
}

#' @exportS3Method base::print
print.Test <- function(x, digits = getOption("digits"), ...) {
  fmt <- function(v) format(v, digits = max(1L, digits - 2L))
  
  cat("\n")
  cat(strwrap(x$name, prefix = "\t"), sep = "\n")
  cat("\n")
  
  out <- paste("statistic =", fmt(x$statistic))
  if (!is.null(x$df)) out <- c(out, paste("df =", fmt(x$df)))
  if (!is.null(x$p_value)) {
    fp <- format.pval(x$p_value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (startsWith(fp, "<")) fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  
  lev <- paste0(format(100 * x$level), "%")
  if (!is.null(x$critical)) {
    cat(sprintf("critical value (%s confidence level) = %s\n", lev, fmt(x$critical)))
  }
  if (!is.null(x$p_value) || !is.null(x$critical)) {
    cat(sprintf("null hypothesis rejected at the %s confidence level: %s\n",
                lev, if (test_reject(x)) "yes" else "no"))
  }
  
  if (length(x$extra)) {
    cat("further statistics:\n")
    print(unlist(x$extra), digits = digits, ...)
  }
  if (nzchar(x$source)) cat("source:", x$source, "\n")
  cat("\n")
  invisible(x)
}



#' Rayleigh test of uniformity against a unimodal alternative in \eqn{S^3}
#' 
#' @details The statistic is \deqn{S = p n \bar{R}^2} for \eqn{n} unit vectors in \eqn{p} dimensions,
#' asymptotically chi-square with \eqn{p} degrees of freedom; the p-value is that of
#' the modified statistic \deqn{S^* = (1 - \frac{1}{2n}) S + \frac{S^2}{2n(p + 2)}} which
#' has error \eqn{O(1/n^2)} (Mardia and Jupp, 2000; eq. 10.4.4-10.4.6, eq. 6.3.3 for the circle; Fisher et al., 1987,
#' sec. 5.3.1(i)). The test is not consistent against axial alternatives.
#'                     
#' @param x Object of class `spherical`
#' @inheritParams bingham_ellipse
#' 
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @references 
#' Mardia, K.V. and Jupp, P.E. (2000) Directional Statistics. John Wiley & Sons, Chichester.
#' 
#' Fisher, N.I., Lewis, T. and Embleton, B.J.J. (1987) Statistical Analysis of Spherical Data. Cambridge University Press, Cambridge.
#' 
#' @export
#' @return object of class `"Test"`
#' @family tests
#' 
#' @examples
#' rayleigh_test3(example_lines)
rayleigh_test3 <- function(x, level = 0.95){
  xv <- Vec3(x)
  d <- dim(xv)
  n <- d[1]
  p <- d[2]
  
  R     <- sqrt(sum(colSums(x)^2))
  S     <- p * R^2 / n
  S_mod <- (1 - 1 / (2 * n)) * S + S^2 / (2 * n * (p + 2))
  
  Test(
      name      = "Rayleigh Test",
      statistic = S,
      p_value   = stats::pchisq(S_mod, df = p, lower.tail = FALSE),
      df        = p,
      critical  = stats::qchisq(level, df = p),
      level     = level,
      source    = "Mardia & Jupp (2000), eq. 10.4.4-10.4.6",
      extra     = list(R = R, Rbar = R / n, modified = S_mod)
    )
}

#' Bingham test of uniformity of axes on the sphere
#' 
#' @details The statistic is \deqn{S = 15 n / 2 (\text{tr}(T^2) - 1/3)} with the orientation 
#' tensor \eqn{T} (scaled to trace 1), asymptotically chi-square with 5 degrees of 
#' freedom (Mardia and Jupp, 2000; eq. 10.7.1-10.7.2), which gives the p-value. 
#' 
#' The modified 
#' statistic \eqn{S^*} of MJ eq. 10.7.3, a polynomial correction meant for \eqn{S}
#' near the critical values, is given in `"extra"`. 
#' 
#' \eqn{S} equals Bingham's (1974, eq. 5.9) 
#' \eqn{X_U^2} and \eqn{n I} with Lisle's (1985) intensity 
#' \deqn{I = 15/2 \sum((\lambda_i - 1/3)^2)}. 
#' 
#' Large \eqn{S} indicates a preferred orientation 
#' (cluster or girdle); a very small \eqn{S} (lower tail) indicates data more regular than random.
#'                     
#' @param D 3x3 orientation tensor (see [ortensor()])
#' @param n integer. Number of vectors
#' @inheritParams rayleigh_test3
#' 
#' @export
#' @return object of class `"Test"`
#' @family tests
#' 
#' @references Mardia, K.V. and Jupp, P.E. (2000) Directional Statistics. John Wiley & Sons, Chichester.
#' 
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @examples
#' t <- example_lines
#' bingham_test(ortensor(t), n = nrow(t))
bingham_test <- function(D, n, level = 0.95) {
  # Bingham test of uniformity of axes on the sphere
  # (Mardia & Jupp 2000, eq. 10.7.1-10.7.3; see the Python docstring).
  S <- 7.5 * n * (sum(D^2) - 1 / 3)
  
  # MJ eq. 10.7.3 for p = 3
  b0 <- 31 / 42
  b1 <- -41 / 294
  b2 <- 20 / 5292
  S_mod <- S * (1 - (b0 + b1 * S + b2 * S^2) / n)
  
  Test(
      name      = "Bingham Test",
      statistic = S,
      p_value   = stats::pchisq(S, df = 5, lower.tail = FALSE),
      df        = 5L,
      critical  = stats::qchisq(level, df = 5),
      level     = level,
      source    = "Mardia & Jupp (2000), eq. 10.7.1-10.7.3",
      extra     = list(modified = S_mod)
    )
}

#' Test of Mean Orientations
#'
#' Test against the null-hypothesis that the samples are drawn from the same Fisher population.
#'
#' @param x,y  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
#'  rows are the observations and the columns are the coordinates.
#' @inheritParams rayleigh_test3
#' @inheritParams sph_mean
#' 
#' @details Given two vector sets \eqn{x} and \eqn{y}, then \eqn{n_x} and 
#' \eqn{n_y} are the number of vectors in each set, \eqn{R_x} and \eqn{R_y} are 
#' the mean resultant lengths of the two vector sets and \eqn{R} is the mean 
#' resultant length of the pooled vector, then 
#' \deqn{F = \frac{(n_x + n_y - 2)(R_x + R_y - R)}{n_x + n_y - R_x - R_y}} 
#' is asymptotically F-distributed with 2 and \eqn{2(n_x + n_y - 2)} degrees of 
#' freedom (Fisher, Lewis & Embleton, 1987, sec. 7.2.2). 
#' Large \eqn{F} indicates that the two samples have different mean orientations.
#' 
#' @references Fisher, N.I., Lewis, T. and Embleton, B.J.J. (1987) Statistical Analysis of Spherical Data. Cambridge University Press, Cambridge.
#'
#' @export
#' @return object of class `"Test"`
#' @family tests
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
fisher_ftest <- function(x, y, level = 0.95, na.rm = TRUE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x) | is.Ray(x))
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
  
  #crit <- stats::qf(p = level, df1 = df1, df2 = df2, lower.tail = FALSE)
  # if (stat > crit) {
  #   message("Reject null-hypothesis")
  # } else {
  #   message("Do not reject null-hypothesis")
  # }
  # c("F stat" = stat, "p-value" = crit)
  
  Test(
    name      = "Fisher F-test for a common mean direction",
    statistic = stat,
    p_value   = stats::pf(stat, df1 = df1, df2 = df2, lower.tail = FALSE),
    df        = c(df1 = df1, df2 = df2),
    critical  = stats::qf(level, df1 = df1, df2 = df2),
    level     = level,
    source    = "Fisher, Lewis & Embleton (1987), sec. 7.2.2",
    extra     = list(Rx = Rx, Ry = Ry, R = R)
  )
}



## ---- helpers --------------------------------------------------------------

.cross3 <- function(a, b) {
  c(a[2L] * b[3L] - a[3L] * b[2L],
    a[3L] * b[1L] - a[1L] * b[3L],
    a[1L] * b[2L] - a[2L] * b[1L])
}

.tangent_basis <- function(mu) {
  # orthonormal e1, e2 spanning the tangent plane at unit vector mu
  ref <- diag(3)[which.min(abs(mu)), ]
  e1 <- .cross3(ref, mu)
  e1 <- e1 / sqrt(sum(e1^2))
  list(e1 = e1, e2 = .cross3(mu, e1))
}

.uniform_scores <- function(phi) {
  # sorted longitudes mapped to [0, 1)
  sort((phi %% (2 * pi)) / (2 * pi))
}

.fisher_kappa_mle <- function(Rbar) {
  # solve A3(kappa) = coth(kappa) - 1/kappa = Rbar
  if (Rbar < 1e-12) return(0)
  if (Rbar >= 1) return(Inf)
  A3 <- function(k) if (k < 1e-4) k / 3 - k^3 / 45 else 1 / tanh(k) - 1 / k
  # k/3 > A3(k) > 1 - 1/k  brackets the root in [3 Rbar, 1 / (1 - Rbar)]
  uniroot(function(k) A3(k) - Rbar,
          lower = 3 * Rbar, upper = 1 / (1 - Rbar),
          tol = .Machine$double.eps^0.75)$root
}

## ---- Fisher goodness of fit ------------------------------------------------

# upper 10, 5 and 1 % points of M_E and M_U (FLE, with eq. 5.15-5.16)
.FISHER_GOF_POINTS <- data.frame(
  level = c(0.90, 0.95, 0.99),
  M_E   = c(0.990, 1.094, 1.308),
  M_U   = c(1.138, 1.207, 1.347)
)

#' Goodness-of-fit tests of the Fisher distribution
#' 
#' @param x object of class `"Spherical"`
#' @inheritParams rayleigh_test3
#' 
#' @details
#'  The data are rotated to their mean direction. Under the Fisher model the
#'  colatitudes \eqn{\theta'} give \eqn{1 - \cos\theta'} approximately exponential
#'  with rate \eqn{k = (n - 1) / \sum(1 - \cos\theta')} and the longitudes are
#'  uniform. 
#'  
#'  The colatitude test \deqn{M_E = (D_n - \frac{0.2}{n})(\sqrt{n} + 0.26 + 
#'  \frac{\sqrt{n}}{2})} uses the Kolmogorov statistic \eqn{D_n} against the
#'  exponential distribution, while the longitude test \deqn{M_U = V_n(\sqrt{n} - 0.467 +
#'  \frac{1.623}{\sqrt{n}})} uses the Kuiper statistic \eqn{V_n} against the uniform 
#'  distribution (Fisher et al., 1987; eq. 5.15-5.16). 
#'  
#'  Critical values are known for the
#'  levels 0.90, 0.95 and 0.99.
#'  
#' @note The constants are FLE's, which are the original ones of Stephens
#'  (1974) and Fisher & Best (1984). Mardia & Jupp (2000, eq. 12.3.1, 
#'  12.3.4) misquote them.
#'  
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @references Fisher, N.I., Lewis, T. and Embleton, B.J.J. (1987) Statistical Analysis of Spherical Data. Cambridge University Press, Cambridge.
#' 
#' @export
#' @return object of class `"Test"`
#' @family tests
#' 
#' @examples
#' set.seed(20250411)
#' x <- rvmf(100)
#' fisher_goodness_of_fit(x)
fisher_goodness_of_fit <- function(x, level = 0.95) {
  X <- Vec3(x)
  row <- which(abs(.FISHER_GOF_POINTS$level - level) < 1e-8)
  if (length(row) != 1L) stop("`level` must be 0.90, 0.95 or 0.99")
  n <- nrow(X)
  
  mu <- colSums(X)
  mu <- mu / sqrt(sum(mu^2))
  B  <- .tangent_basis(mu)
  i  <- seq_len(n)
  
  # colatitudes: Kolmogorov D_n against the exponential
  x <- sort(1 - pmin(pmax(drop(X %*% mu), -1), 1))
  k <- (n - 1) / sum(x)
  F <- 1 - exp(-k * x)
  D <- max(i / n - F, F - (i - 1) / n)
  m_e <- (D - 0.2 / n) * (sqrt(n) + 0.26 + 0.5 / sqrt(n))
  
  # longitudes: Kuiper V_n against the uniform
  phi <- atan2(drop(X %*% B$e2), drop(X %*% B$e1))
  u <- .uniform_scores(phi)
  V <- max(i / n - u) + max(u - (i - 1) / n)
  m_u <- V * (sqrt(n) - 0.467 + 1.623 / sqrt(n))
  
  source <- "Fisher, Lewis & Embleton (1987), eq. 5.15-5.16"
  list(
    colatitude = Test(
      "Fisher colatitude M_E", m_e,
      critical = .FISHER_GOF_POINTS$M_E[row], level = level,
      source = source, extra = list(D = D, k = k)
    ),
    longitude = Test(
      "Fisher longitude M_U", m_u,
      critical = .FISHER_GOF_POINTS$M_U[row], level = level,
      source = source, extra = list(V = V)
    )
  )
}

## ---- Fisher vs Kent score test ---------------------------------------------

#' Score test of the Fisher against the Kent distribution
#' 
#' @details \deqn{T = n \frac{(\kappa/2)^2 I_{1/2}(\kappa)}{I_{5/2}(\kappa)} (r_2^2)} is asymptotically 
#' \eqn{\chi^2} with 2 degrees of freedom (Kent, 1982, eq. 5.5), with the 
#' maximum-likelihood \eqn{\kappa} of the Fisher distribution (Kent's Theorem 
#' 5.1) and \eqn{r_2} the difference of the two eigenvalues of the second 
#' moments of the data in the tangent plane of the mean direction. 
#' 
#' Large \eqn{T}
#' indicates elliptical (Kent) rather than circular (Fisher) contours.
#' 
#' @note
#' FLE (sec. 5.4.6) computes the statistic with the concentration of the 
#' fitted Kent distribution, which gives a larger value; this follows 
#' Kent's theorem.
#' 
#' @references 
#' Kent, J.T. (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society B 44, 71-80.
#' 
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @inherit fisher_goodness_of_fit 
#' 
#' @export
#' @return object of class `"Test"`
#' @family tests
#' 
#' @examples
#' #' set.seed(20250411)
#' x <- rvmf(100)
#' fisher_kent_test(x)
fisher_kent_test <- function(x, level = 0.95) {
  X <- Vec3(x)
  n <- nrow(X)
  resultant <- colSums(X)
  R <- sqrt(sum(resultant^2))
  Rbar <- R / n
  
  B <- .tangent_basis(resultant / R)
  P <- X %*% cbind(B$e1, B$e2)
  lam <- eigen(crossprod(P) / n, symmetric = TRUE, only.values = TRUE)$values
  r2 <- lam[1L] - lam[2L]                       # eigen() sorts descending
  
  kappa <- .fisher_kappa_mle(Rbar)
  T <- n * (kappa / 2)^2 *
    besselI(kappa, 0.5, expon.scaled = TRUE) /
    besselI(kappa, 2.5, expon.scaled = TRUE) * r2^2
  
  Test(
    name      = "Fisher vs Kent",
    statistic = T,
    p_value   = stats::pchisq(T, df = 2, lower.tail = FALSE),
    df        = 2L,
    critical  = stats::qchisq(level, df = 2),
    level     = level,
    source    = "Kent (1982), eq. 5.5",
    extra     = list(kappa = kappa, r2 = r2)
  )
}

## ---- Watson symmetry tests -------------------------------------------------


#' Tests of rotational symmetry of a Bingham distribution
#' 
#' @details Whether a Bingham distribution is a Watson distribution: 
#' \deqn{w_b = n \sum((k_i - \text{mean}(k)) (\lambda_i - \text{mean}(\lambda)))} 
#' over the two smallest eigenvalues tests the bipolar form, \eqn{w_g} over the 
#' two largest the girdle form, each asymptotically chi-square with 2 degrees of 
#' freedom (Mardia and Jupp, 2000; eq. 10.7.8-10.7.18). \eqn{k_i} are the maximum-likelihood 
#' Bingham concentrations \eqn{(0, \kappa_\text{major}, \kappa_\text{minor})} and
#' \eqn{\lambda_i} the eigenvalues of the orientation tensor, both in descending order. 
#' Mardia & Jupp's calcite example gives \eqn{w_g = 1.33} (girdle accepted) and \eqn{w_b = 9.97}.
#'  
#' @param lambda numeric.
#' @param k Bingham concentrations \eqn{(0, \kappa_\text{major}, \kappa_\text{minor})} 
#' @inheritParams bingham_test
#' 
#' @references Mardia K.V. and Jupp P.E. (2000) Directional Statistics. John Wiley & Sons, Chichester.
#' 
#' @export
#' @return object of class `"Test"`
#' @family tests
#' 
#' @examples
#' r <- rbingham_eig(100, c(100, 1, 0), "Line")
#' e <- ot_eigen(r)$values
#' k <- c(0, bingham_statistics(r)$k)
#' watson_symmetry_test(e, k, n = nrow(r))
watson_symmetry_test <- function(lambda, k, n, level = 0.95) {
  t <- as.numeric(lambda)
  k <- as.numeric(k)
  idx <- list(bipolar = 2:3, girdle = 1:2)      # 1-based
  lapply(setNames(names(idx), names(idx)), function(nm) {
    j <- idx[[nm]]
    w <- n * sum((k[j] - mean(k[j])) * (t[j] - mean(t[j])))
    Test(
      name      = sprintf("Watson %s symmetry", nm),
      statistic = w,
      p_value   = stats::pchisq(w, df = 2, lower.tail = FALSE),
      df        = 2L,
      critical  = stats::qchisq(level, df = 2),
      level     = level,
      source    = "Mardia & Jupp (2000), eq. 10.7.8-10.7.18"
    )
  })
}


# Comparing samples ------------------------------------------------------------

samples <- list(
  x1 = Line(rvmf()), x2 = Line(rvmf())
)

.resultants <- function(samples) {
  samples <- lapply(samples, Vec3)
  
  # Sizes, resultant lengths and resultant vectors of a list of unit-vector
  # matrices (each n_i x p).
  if (!is.list(samples) || length(samples) < 2L) {
    stop("At least two samples are needed")
  }
  samples <- lapply(samples, function(x) {
    d <- dim(x)
    matrix(as.numeric(x), d[1L], d[2L])
  })
  p <- unique(vapply(samples, ncol, integer(1L)))
  if (length(p) != 1L) stop("The samples must have the same dimension")
  
  n <- vapply(samples, nrow, integer(1L))
  S <- do.call(rbind, lapply(samples, colSums))    # k x p, one resultant per row
  list(n = n, R = sqrt(rowSums(S^2)), S = S, p = p)
}

.f_result <- function(name, stat, df, level, source) {
  Test(
    name      = name,
    statistic = stat,
    p_value   = pf(stat, df[1L], df[2L], lower.tail = FALSE),
    df        = c(df1 = df[1L], df2 = df[2L]),
    critical  = qf(level, df[1L], df[2L]),
    level     = level,
    source    = source
  )
}

.chi2_result <- function(name, stat, df, level, source) {
  Test(
    name      = name,
    statistic = stat,
    p_value   = pchisq(stat, df, lower.tail = FALSE),
    df        = df,
    critical  = qchisq(level, df),
    level     = level,
    source    = source
  )
}

.vonmises_kappa_mle <- function(Rbar) {
  # solve A1(kappa) = I1(kappa) / I0(kappa) = Rbar
  if (Rbar < 1e-12) return(0)
  if (Rbar >= 1) return(Inf)
  A1 <- function(k) besselI(k, 1, expon.scaled = TRUE) /
    besselI(k, 0, expon.scaled = TRUE)
  uniroot(function(k) A1(k) - Rbar, lower = 0, upper = 1 / (1 - Rbar),
          extendInt = "upX", tol = .Machine$double.eps^0.75)$root
}

### ---- common mean test -----------------------------------------------------

common_mean_test <- function(samples, method = NULL, level = 0.95) {
  res <- .resultants(samples)
  n_i <- res$n
  R_i <- res$R
  S_i <- res$S
  p   <- res$p
  if (!p %in% 2:3) stop("The samples must be 2D (circle) or 3D (sphere)")
  
  n    <- sum(n_i)
  q    <- length(n_i)
  sumR <- sum(R_i)
  R    <- sqrt(sum(colSums(S_i)^2))
  
  sphere <- c("watson", "mcfadden_lowes", "watson_debiche")
  circle <- c("watson_williams", "uniform_scores", "watson_u2")
  allowed <- if (p == 3L) sphere else circle
  if (is.null(method)) method <- allowed[1L]
  if (!method %in% allowed) {
    stop(sprintf("Unknown method '%s' for %dD data, expected %s",
                 method, p, paste(sort(allowed), collapse = ", ")))
  }
  
  if (method %in% c("watson", "watson_williams")) {
    F <- (n - q) * (sumR - R) / ((q - 1) * (n - sumR))
    if (p == 3L) {
      df <- c(2 * (q - 1), 2 * (n - q))
      source <- "Watson (1956), eq. 2.32; Mardia & Jupp (2000), eq. 10.6.4"
    } else {
      kappa <- .vonmises_kappa_mle(sumR / n)
      F <- F * (if (kappa > 0) 1 + 3 / (8 * kappa) else Inf)
      df <- c(q - 1, n - q)
      source <- "Mardia & Jupp (2000), eq. 7.4.4-7.4.6"
    }
    return(.f_result(paste(method, "common mean"), F, df, level, source))
  }
  
  if (method == "mcfadden_lowes") {
    g <- (n - q) / (q - 1) * (sumR - R^2 / sumR) / (2 * (n - sumR))
    return(.f_result("McFadden-Lowes common mean", g,
                     c(2 * (q - 1), 2 * (n - q)), level,
                     "McFadden & Lowes (1981), eq. 25, 41"))
  }
  
  if (method == "watson_debiche") {
    k <- (n_i - 1) / (n_i - R_i)
    V <- 2 * (sum(k * R_i) - sqrt(sum(colSums(k * S_i)^2)))
    return(.chi2_result("Watson-Debiche common mean", V, 2 * (q - 1), level,
                        "Watson & Debiche (1992), eq. 22, 25-26"))
  }
  
  theta <- lapply(samples, function(x) {
    x <- matrix(as.numeric(x), ncol = 2L)
    atan2(x[, 2L], x[, 1L])
  })
  if (method == "uniform_scores") {
    return(.uniform_scores_test(theta, level))
  }
  .watson_u2_two_sample(theta, level)
}


### ---- circular two-sample tests --------------------------------------------

.uniform_scores_test <- function(theta, level) {
  pooled <- unlist(theta, use.names = FALSE) %% (2 * pi)
  n <- length(pooled)
  groups <- rep(seq_along(theta), lengths(theta))
  # ties within a sample get consecutive ranks, ties between samples the
  # average rank
  ranks <- rank(pooled, ties.method = "first")
  for (value in unique(pooled[duplicated(pooled)])) {
    tied <- pooled == value
    if (length(unique(groups[tied])) > 1L) ranks[tied] <- mean(ranks[tied])
  }
  beta <- 2 * pi * ranks / n
  CS  <- rowsum(cbind(cos(beta), sin(beta)), groups)   # q x 2: C_i, S_i
  n_i <- tabulate(groups)
  W <- if (length(theta) == 2L) {
    2 * (n - 1) * sum(CS[1L, ]^2) / (n_i[1L] * n_i[2L])
  } else {
    2 * sum(rowSums(CS^2) / n_i)
  }
  .chi2_result("uniform scores", W, 2 * (length(theta) - 1), level,
               "Mardia & Jupp (2000), eq. 8.3.3-8.3.4, 8.4.3")
}

.watson_u2_sf <- function(u, m_max = 1000L) {
  # asymptotic P(U^2 > u) = 2 sum_{m>=1} (-1)^(m-1) exp(-2 m^2 pi^2 u)
  if (u <= 0) return(1)
  m <- seq_len(m_max)
  min(max(2 * sum((-1)^(m - 1) * exp(-2 * m^2 * pi^2 * u)), 0), 1)
}

.watson_u2_isf <- function(alpha) {
  uniroot(function(u) .watson_u2_sf(u) - alpha, lower = 1e-3, upper = 5,
          tol = .Machine$double.eps^0.75)$root
}

.watson_u2_two_sample <- function(theta, level) {
  if (length(theta) != 2L) stop("Watson's two-sample U2 compares two samples")
  n1 <- length(theta[[1L]])
  n2 <- length(theta[[2L]])
  n  <- n1 + n2
  pooled <- unlist(theta, use.names = FALSE) %% (2 * pi)
  first <- rep(c(1, 0), c(n1, n2))[order(pooled)]     # order() is stable
  a <- cumsum(first)
  d <- (seq_len(n) - a) / n2 - a / n1
  u2 <- n1 * n2 / n^2 * sum((d - mean(d))^2)
  Test(
    name      = "Watson two-sample U2",
    statistic = u2,
    p_value   = .watson_u2_sf(u2),
    critical  = .watson_u2_isf(1 - level),
    level     = level,
    source    = "Mardia & Jupp (2000), eq. 8.3.7"
  )
}

### ---- equal concentration --------------------------------------------------

equal_kappa_test <- function(samples, method = NULL, level = 0.95) {
  res <- .resultants(samples)
  n_i <- res$n
  q   <- length(n_i)
  nu  <- (res$p - 1) * (n_i - 1)
  dev <- n_i - res$R
  if (is.null(method)) method <- if (q == 2L) "f" else "bartlett"
  
  if (method == "f") {
    if (q != 2L) stop("The F test compares two samples")
    F <- (dev[1L] / nu[1L]) / (dev[2L] / nu[2L])
    p_value <- 2 * min(pf(F, nu[1L], nu[2L]),
                       pf(F, nu[1L], nu[2L], lower.tail = FALSE))
    return(Test(
      name      = "equal concentration F",
      statistic = F,
      p_value   = min(p_value, 1),
      df        = c(df1 = nu[1L], df2 = nu[2L]),
      level     = level,
      source    = "Mardia & Jupp (2000), eq. 10.5.3"
    ))
  }
  if (method == "bartlett") {
    total <- sum(nu)
    B <- total * log(sum(dev) / total) - sum(nu * log(dev / nu))
    B <- B / (1 + (sum(1 / nu) - 1 / total) / (3 * (q - 1)))
    return(.chi2_result("equal concentration Bartlett", B, q - 1, level,
                        "Mardia & Jupp (2000), eq. 10.6.21"))
  }
  stop(sprintf("Unknown method '%s', expected f or bartlett", method))
}

### ---- axial common mean ----------------------------------------------------

axial_common_mean_test <- function(samples, level = 0.95) {
  if (!is.list(samples) || length(samples) < 2L) {
    stop("At least two samples are needed")
  }
  samples <- lapply(samples, function(x) matrix(as.numeric(x), ncol = 3L))
  taus <- vapply(samples, function(x) {
    eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values[1L]
  }, numeric(1L))
  
  X <- do.call(rbind, samples)
  e <- eigen(crossprod(X), symmetric = TRUE)
  tau_bar <- e$values[1L] / nrow(X)
  G <- mean(drop(X %*% e$vectors[, 1L])^4)
  g <- (tau_bar - G) / (3 * tau_bar - 1)
  N <- (sum(taus) - e$values[1L]) / g
  
  .chi2_result("axial common mean", N, 2 * (length(samples) - 1), level,
               "Fisher, Lewis & Embleton (1987), eq. 7.72-7.76")
}

### ---- frame uniformity -----------------------------------------------------

frame_uniformity_test <- function(frames, level = 0.95) {
  # frames: n x r x p array, frames[k, j, ] is axis j of frame k
  d <- dim(frames)
  if (length(d) != 3L) stop("`frames` must be an n x r x p array")
  U <- array(as.numeric(frames), d)
  n <- d[1L]; r <- d[2L]; p <- d[3L]
  I <- diag(p) / p
  
  T <- lapply(seq_len(r), function(j) {
    crossprod(matrix(U[, j, ], n, p)) / n - I
  })
  total <- sum(vapply(T, function(Tj) sum(Tj^2), numeric(1L)))
  if (r < p) total <- total + sum(Reduce(`+`, T)^2) / (p - r)
  
  S  <- n * (p - 1) * (p + 2) / 2 * total
  df <- min(r, p - 1) * (p - 1) * (p + 2) / 2
  d1 <- sum(vapply(T, function(Tj) 1 - sum((Tj + I)^2), numeric(1L)))
  
  Test(
    name      = "Arnold-Jupp frame uniformity",
    statistic = S,
    p_value   = pchisq(S, df, lower.tail = FALSE),
    df        = df,
    critical  = qchisq(level, df),
    level     = level,
    source    = "Arnold & Jupp (2013), eq. 5",
    extra     = list(dispersion = d1)
  )
}