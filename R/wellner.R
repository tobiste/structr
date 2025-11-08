### WELLNER (1979) TWO-SAMPLE TEST ###

#' Wellner's Two-Sample Test
#'
#' Wellner's (1979, Example 1a) Rayleigh-style T-statistic, which 
#' quantifies the dissimilarity between two sets of vectors. The statistic 
#' increases with the degree of difference between the datasets.
#'
#' `wellner_inference()` performs a permutation-based inference using 
#' Wellner's T-statistic to assess whether the two datasets are drawn from 
#' the same population.
#'
#' @inheritParams reject
#' @param n_perm integer. Number of permutations.
#' 
#' @return 
#' `wellner()` computes Wellner's T-statistic, a non-negative measure of 
#' dissimilarity between two datasets. The value is zero when the datasets 
#' are identical.
#'
#' `wellner_inference()` estimates the fraction of permutation tests in which 
#' the computed T-statistic exceeds the observed T for the original data 
#' (a number between 0 and 1). This value can be interpreted as a p-value 
#' for the null hypothesis that the two populations are identical 
#' (not merely that their means coincide). 
#' Thus, smaller p-values indicate stronger evidence that the two populations 
#' differ in a meaningful way.
#' 
#' @references  Jon A. Wellner. "Permutation Tests for Directional Data." 
#' Ann. Statist. 7(5) 929-943, September, 1979. \doi{10.1214/aos/1176344779}
#' 
#' @name wellner
#' 
#' @examples
#' test <- rvmf(100)
#' wellner(test, Line(120, 50))
#' wellner_inference(test, Line(120, 50))
NULL


#' @rdname wellner
#' @export
wellner <- function(x, y) UseMethod("wellner")
  
#' @rdname wellner
#' @export
wellner.Ray <- function(x, y) {
  xs <- vec_list(x)
  ys <- vec_list(y)
  
  # The number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(xs)
  n <- length(ys)
  p <- length(xs[[1]]) - 1
  # Compute Wellner's statistic.
  diff <- arithmeticMean(xs) - arithmeticMean(ys)
  dot(diff, diff) * (p + 1) * m * n / (m + n)
}

#' @rdname wellner
#' @export
wellner.Line <- function(x, y) {
  xs <- vec_list(x)
  ys <- vec_list(y)
  
  # The number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(xs)
  n <- length(ys)
  p <- length(xs[[1]]) - 1
  # Compute m Tx and n Ty.
  f <- function(v) {
    outer(v, v)
  }
  mTx <- Reduce("+", lapply(xs, f))
  nTy <- Reduce("+", lapply(ys, f))
  # Compute Wellner's statistic.
  txMinusTy <- mTx / m - nTy / n
  tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}

#' @rdname wellner
#' @export
wellner.Vec3 <- function(x, y) wellner.Line(x, y)

#' @rdname wellner
#' @export
wellner.Plane <- function(x, y) wellner.Line(x, y)


# Also computes Wellner's T-statistic, with inputs presented in a manner suitable for fast permutation testing.
wellner_shortcut_ray <- function(xys, choices, rxPlusRy) {
  # Prepare the number of xs and ys and the dimension of the hypersphere.
  m <- length(choices)
  n <- length(xys) - m
  p <- length(xys[[1]]) - 1
  # Compute Rx and Ry.
  rx <- Reduce("+", lapply(choices, function(i) xys[[i]]))
  ry <- rxPlusRy - rx
  # Compute Wellner's statistic.
  diff <- rx / m - ry / n
  dot(diff, diff) * (p + 1) * m * n / (m + n)
}

#' Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).
#'
#' Assumes large sample sizes, specifically that choose(nx + ny, nx) >> numPerms. For small sample sizes, see rayWellnerExactInference.
#' @param xs A list of rays.
#' @param ys A list of rays.
#' @param numPerms A real number (positive integer). The number of permutations, say 1,000 or 10,000.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
Wellner_inference_ray <- function(xs, ys, numPerms) {
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  rxPlusRy <- Reduce("+", xys)
  # Compute Wellner's statistic for the actual data.
  t <- wellner_shortcut_ray(xys, 1:m, rxPlusRy)
  # Compute Wellner's statistic for permuted data, numPerms times.
  ts <- replicate(
    numPerms, wellner_shortcut_ray(xys, sample.int(m + n, m), rxPlusRy)
  )
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function(u) {
    u > t
  }, ts)
  length(greaterThan) / numPerms
}

#' Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).
#'
#' Assumes small sample sizes. Deterministically generates all choose(nx + ny, nx) reassignments of the data to the two groups. For large sample sizes, see rayWellnerInference.
#' @param xs A list of rays.
#' @param ys A list of rays.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
wellner_exact_inference_ray <- function(xs, ys) {
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  rxPlusRy <- Reduce("+", xys)
  # Compute Wellner's statistic for the actual data.
  t <- wellner_shortcut_ray(xys, 1:m, rxPlusRy)
  # Compute Wellner's statistic for all possible data permutations.
  rows <- combs(1:(m + n), m)
  ts <- sapply(
    1:nrow(rows),
    function(i) wellner_shortcut_ray(xys, rows[i, ], rxPlusRy)
  )
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function(u) {
    u > t
  }, ts)
  length(greaterThan) / choose(m + n, m)
}



# Also computes Wellner's T-statistic, with inputs presented in a manner suitable for fast permutation testing.
wellner_shortcut_line <- function(xys, choices, mTxPlusnTy) {
  # The number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(choices)
  n <- length(xys) - m
  p <- length(xys[[1]]) - 1
  # Compute m Tx and n Ty.
  mTx <- Reduce("+", lapply(choices, function(i) {
    outer(xys[[i]], xys[[i]])
  }))
  nTy <- mTxPlusnTy - mTx
  # Compute Wellner's statistic.
  txMinusTy <- mTx / m - nTy / n
  tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}


#' @rdname wellner
#' @export
wellner_inference <- function(x, y, n_perm) UseMethod("wellner_inference")

#' @rdname wellner
#' @export
wellner_inference.Vec3 <- function(x, y, n_perm = 1000) wellner_inference.Line(x, y, n_perm = 1000)

#' @rdname wellner
#' @export
wellner_inference.Line <- function(x, y, n_perm = 1000) {
  xs <- vec_list(x)
  ys <- vec_list(y)
  
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  mTxPlusnTy <- Reduce("+", lapply(xys, function(v) {
    outer(v, v)
  }))
  # Compute Wellner's statistic for the actual data.
  t <- wellner_shortcut_line(xys, 1:m, mTxPlusnTy)
  # Compute Wellner's statistic for permuted data, n_perm times.
  ts <- replicate(
    n_perm, wellner_shortcut_line(xys, sample.int(m + n, m), mTxPlusnTy)
  )
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function(u) {
    u > t
  }, ts)
  length(greaterThan) / n_perm
}

#' @rdname wellner
#' @export
wellner_inference.Ray <- function(x, y, n_perm = 1000) {
  xs <- vec_list(x)
  ys <- vec_list(y)
  
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  rxPlusRy <- Reduce("+", xys)
  # Compute Wellner's statistic for the actual data.
  t <- wellner_shortcut_ray(xys, 1:m, rxPlusRy)
  # Compute Wellner's statistic for permuted data, n_perm times.
  ts <- replicate(
    n_perm, wellner_shortcut_ray(xys, sample.int(m + n, m), rxPlusRy)
  )
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function(u) {
    u > t
  }, ts)
  length(greaterThan) / n_perm
}

#' @rdname wellner
#' @export
wellner_inference.Plane <- function(x, y, n_perm = 1000) wellner_inference.Line(x, y, n_perm = 1000)
