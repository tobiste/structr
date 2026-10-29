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



#' Rayleigh test of uniformity against a unimodal alternative
#' 
#' @details The statistic is \eqn{S = p n \bar{R}^2} for \eqn{n} unit vectors in \eqn{p} dimensions,
#' asymptotically chi-square with \eqn{p} degrees of freedom; the p-value is that of
#' the modified statistic \eqn{S^* = (1 - 1/(2n)) S + S^2 / (2n(p + 2))}, which
#' has error \eqn{O(1/n^2)} (MJ eq. 10.4.4-10.4.6, eq. 6.3.3 for the circle; FLE
#' sec. 5.3.1(i)). The test is not consistent against axial alternatives.
#'                     
#' @param x Object of class `spherical`
#' @inheritParams bingham_ellipse
#' 
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
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
#' freedom (MJ eq. 10.7.1-10.7.2), which gives the p-value. 
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
#' @source modified after the python package `apsg` (v.2.0.6) by O. Lexa
#' 
#' @family tests
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
#' @return object of class `"Test"`
#'
#' @export
#'
#' @importFrom stats qf
#' 
#' @family tests
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
#'  colatitudes \eqn{\theta'} give \eqn{1 - \cos(\theta')} approximately exponential
#'  with rate \eqn{k = (n - 1) / \sum(1 - \cos(\theta'))} and the longitudes are
#'  uniform. 
#'  
#'  The colatitude test \deqn{M_E = (D_n - 0.2/n)(\sqrt(n) + 0.26 + 
#'  0.5/\sqrt(n))} uses the Kolmogorov statistic \eqn{D_n} against the
#'  exponential distribution, while the longitude test \deqn{M_U = V_n(\sqrt(n) - 0.467 +
#'  1.623/\sqrt(n))} uses the Kuiper statistic \eqn{V_n} against the uniform 
#'  distribution (FLE eq. 5.15-5.16). 
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
#' @return object of class `"Test"`
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
#' @details \eqn{T = n (\kappa/2)^2 I_1/2(\kappa) / I_5/2(\kappa) r_2^2} is asymptotically 
#' chi-square with 2 degrees of freedom (Kent, 1982, eq. 5.5), with the 
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
#' @return object of class `"Test"`
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
#' freedom (MJ eq. 10.7.8-10.7.18). \eqn{k_i} are the maximum-likelihood 
#' Bingham concentrations \eqn{(0, \kappa_\text{major}, \kappa_\text{minor})} and
#' \eqn{\lambda_i} the eigenvalues of the orientation tensor, both in descending order. 
#' Mardia & Jupp's calcite example gives \eqn{w_g = 1.33} (girdle accepted) and \eqn{w_b = 9.97}.
#'  
#' @param lambda numeric.
#' @param k Bingham concentrations \eqn{(0, \kappa_\text{major}, \kappa_\text{minor})} 
#' @inheritParams bingham_test
#' 
#' @return object of class `"Test"`
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
