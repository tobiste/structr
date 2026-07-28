# WISSI: Weighted Iterative Sigma-Space Inversion
#
# A fault slip stress inversion method combining:
#   - Michael (1984)       : analytic uncertainty framework
#   - Angelier (1990)      : upsilon magnitude-aware misfit criterion
#   - Yamaji & Sato (2006) : 5-sphere sigma-space, angular stress distance
#   - Hansen (2013)        : robust handling of unknown/uncertain slip sense
#
# References
# ---
#   Michael, A.J. (1984). J. Geophys. Res., 89, 11517-11526.
#   Angelier, J. (1990). J. Geophys. Res., 95, 17357-17383.
#   Mostafa, M.E. (2005). J. Struct. Geol., 27, 930-940.
#   Yamaji, A. & Sato, K. (2006). Geophys. J. Int., 167, 933-942.
#   Pascal, C. (2022). Paleostress Inversion Techniques. Elsevier.
#
# Coordinate convention
# ---
#   Right-handed Cartesian: x = East, y = North, z = Up.
#   normals : N x 3 matrix of unit fault plane normals (upward-pointing).
#   slips   : N x 3 matrix of unit slip direction vectors, parallel to the
#             shear traction (direction of hanging-wall motion).
#
# Opposite-tensor note
# ---
#   The Wallace-Bott condition is satisfied by both y and -y, corresponding
#   to tensors with swapped sigma1/sigma3 and Phi <-> 1-Phi. The sense
#   annealing (Stage 3) resolves this when slip senses are known. When all
#   senses are unknown, both solutions are returned; the geologist selects
#   based on independent evidence (fault kinematics, regional context).
# =============================================================================


# ===
# INTERNAL HELPERS  (prefixed with . — not intended for direct use)
# ===

# Build normalised reduced stress tensor from Phi and principal axes
# R_mat: 3x3, columns = (sigma1, sigma2, sigma3) unit direction vectors
.build_sigma_from_phi <- function(Phi, R_mat) {
  Rsc <- sqrt(3 * Phi^2 - 3 * Phi + 3)
  R_mat %*% diag(c(2 - Phi, 2 * Phi - 1, -Phi - 1) / Rsc) %*% t(R_mat)
}

# Shear traction vector on fault plane i
.shear_tau <- function(sigma, n) {
  t <- sigma %*% n
  drop(t - sum(t * n) * n)
}

# Shear traction magnitude
.tau_mag <- function(sigma, n) sqrt(sum(.shear_tau(sigma, n)^2))


# Build the 5x5 M matrix from weighted eps' vectors
.build_M5 <- function(normals, slips, weights) {
  N <- nrow(normals)
  M5 <- matrix(0, 5, 5)
  for (i in seq_len(N)) {
    ep5 <- drop(crossprod(.B, .eps6_prime(slips[i, ], normals[i, ])))
    M5 <- M5 + weights[i] * tcrossprod(ep5)
  }
  M5
}

# Solve the 5x5 eigenproblem: return smallest-eigenvalue eigenvector as y6
.solve_M5 <- function(M5) {
  eig <- eigen(M5, symmetric = TRUE) # eigenvalues decreasing
  x5 <- eig$vectors[, 5L]
  y <- drop(.B %*% x5)
  list(
    y = y / sqrt(sum(y^2)),
    eigvals = eig$values,
    vectors = eig$vectors
  )
}



# Per-fault angular misfit alpha (unsigned 0-90 deg, line metric)
.alpha_deg <- function(sigma, normals, slips) {
  vapply(seq_len(nrow(normals)), function(i) {
    tau <- .shear_tau(sigma, normals[i, ])
    tn <- sqrt(sum(tau^2))
    if (tn < 1e-12) {
      return(NA_real_)
    }
    acos(pmax(-1, pmin(1, abs(sum(tau / tn * slips[i, ]))))) * 180 / pi
  }, numeric(1L))
}

# Per-fault signed misfit (0-180 deg); > 90 flags potential sense errors
.alpha_signed_deg <- function(sigma, normals, slips) {
  vapply(seq_len(nrow(normals)), function(i) {
    tau <- .shear_tau(sigma, normals[i, ])
    tn <- sqrt(sum(tau^2))
    if (tn < 1e-12) {
      return(NA_real_)
    }
    acos(pmax(-1, pmin(1, sum(tau / tn * slips[i, ])))) * 180 / pi
  }, numeric(1L))
}

# Extract Phi and principal axes from y6
.decompose_y <- function(y) {
  sigma <- .y6_to_sigma(y)
  eig_s <- eigen(sigma, symmetric = TRUE)
  eigval <- eig_s$values
  eigvec <- eig_s$vectors
  Phi <- (eigval[1] - eigval[2]) / (eigval[1] - eigval[3])
  list(
    sigma = sigma, sigma1 = eigvec[, 1L], sigma2 = eigvec[, 2L],
    sigma3 = eigvec[, 3L], eigenvalues = eigval, Phi = Phi
  )
}



# ===
# WEIGHT UTILITIES
# ===

#' Convert a raw signal vector to weights with NA imputation.
#'
#' Imputation is done before scaling, in the raw signal space.
#' Normalisation is deliberately omitted — done by combine_weights() or
#' internally by wissi().
#'
#' @param signal    Numeric vector (e.g. field ranks 1-5, errors in degrees).
#'                  May contain NAs.
#' @param scale_fn  Function mapping imputed signal -> raw weights.
#'                  Default: 1/x^2 (inverse variance).
#' @param pessimistic If TRUE, impute NAs with the 75th percentile of the
#'                  observed signal (use when missingness is informative).
#' @noRd 
NULL

# signal_to_weights <- function(signal,
#                               scale_fn = function(x) 1 / x^2,
#                               pessimistic = FALSE) {
#   stopifnot(is.numeric(signal), any(!is.na(signal)))
#   fill <- if (pessimistic) {
#     quantile(signal, 0.75, na.rm = TRUE)
#   } else {
#     mean(signal, na.rm = TRUE)
#   }
#   complete <- ifelse(is.na(signal), fill, signal)
#   scale_fn(complete)
# }

# Huber kernel: linear penalty beyond threshold k (default k = 1.5 sigma_mad)
.w_huber <- function(alpha_deg, k = 30) {
  ifelse(alpha_deg <= k, 1, k / pmax(alpha_deg, 1e-6))
}

# Tukey bisquare: hard zeroes beyond threshold (most aggressive outlier rejection)
.w_tukey <- function(alpha_deg, k = 45) {
  ifelse(alpha_deg >= k, 0, (1 - (alpha_deg/k)^2)^2)
}

# Andrews sine: smooth zero beyond threshold
.w_andrews <- function(alpha_deg, k = 60) {
  ifelse(alpha_deg >= k, 0,
         sin(pi * alpha_deg / k) / (pi * alpha_deg / k + 1e-15))
}

# Kernel weight function
.w_robust <- function(robust_kernel, robust_k_deg){
  switch(robust_kernel,
                    huber    = function(a) .w_huber(a,   k = robust_k_deg),
                    tukey    = function(a) .w_tukey(a,   k = robust_k_deg),
                    andrews  = function(a) .w_andrews(a, k = robust_k_deg)
  )
  }

# Spatial median on S5: more robust initial estimate than mean eigenvector
.s5_median <- function(normals, slips, weights, max_iter = 50L) {
  N   <- nrow(normals)
  # Single-fault y-vectors (poles of solution great circles)
  poles <- matrix(0, N, 5)
  for (i in seq_len(N)) {
    ep5 <- drop(crossprod(.B, .eps6_prime(slips[i,], normals[i,])))
    poles[i,] <- ep5 / sqrt(sum(ep5^2))
  }
  # Weiszfeld algorithm on S5
  y_cur <- colMeans(poles * weights)
  y_cur <- y_cur / sqrt(sum(y_cur^2))
  for (iter in seq_len(max_iter)) {
    # Distances from current estimate
    dots  <- pmax(1e-6, abs(poles %*% y_cur))
    theta <- acos(pmin(1, dots))               # angular distances
    # Weiszfeld weights: 1/theta (downweight distant points)
    wz    <- weights / pmax(theta, 1e-4)
    y_new <- colSums(poles * wz)
    y_new <- y_new / sqrt(sum(y_new^2))
    if (acos(pmax(-1,pmin(1,abs(sum(y_cur*y_new))))) * 180/pi < 1e-4) break
    y_cur <- y_new
  }
  y_cur
}

# Flag outliers using MAD-based threshold
.flag_outliers_mad <- function(alpha_deg, k = 2.5) {
  med  <- median(alpha_deg, na.rm = TRUE)
  mad  <- median(abs(alpha_deg - med), na.rm = TRUE)
  # Scale factor 1.4826 makes MAD consistent with normal distribution sigma
  sigma_mad <- 1.4826 * mad
  alpha_deg > (med + k * sigma_mad)
}

# Combine multiple weight vectors multiplicatively and normalise.
#
# Each element of ... is a numeric vector of length N.
# Returns a normalised weight vector (mean = 1).
# combine_weights <- function(...) {
#   wlist <- list(...)
#   if (length(wlist) == 0L) stop("Provide at least one weight vector.")
#   lens <- lengths(wlist)
#   if (length(unique(lens)) != 1L) {
#     stop("All weight vectors must have the same length.")
#   }
#   w <- Reduce(`*`, wlist)
#   if (any(w < 0)) stop("Combined weights must be non-negative.")
#   if (!any(w > 0)) stop("At least one weight must be positive.")
#   w / mean(w)
# }


# ===
# STAGE 1 — Weighted eigenvector initialisation
# ===
# Builds M5 = sum_i omega_i * ep5_i ep5_i^T and returns the
# smallest-eigenvalue eigenvector as the initial stress estimate.
# Equivalent to a weighted Yamaji-Sato inversion.
#
# Returns: list(y, eigvals, M5, eigenvalue_gap)
# ===
.wissi_stage1 <- function(normals, slips, weights) {
  M5 <- .build_M5(normals, slips, weights)
  res <- .solve_M5(M5)
  y <- .resolve_sense(res$y, normals, slips)
  gap <- res$eigvals[4L] - res$eigvals[5L]
  list(
    y = y, eigvals = res$eigvals, M5 = M5,
    eigenvalue_gap = gap
  )
}


# ===
# STAGE 2 — Magnitude reweighting (Angelier/Mostafa in sigma-space)
# ===
# Iteratively replaces global lambda with per-fault mu_i = |tau_i|/lambda_max,
# downweighting mechanically degenerate fault planes. Convergence is measured
# by the angular stress distance Theta between successive iterates — a
# physically meaningful criterion (Theta ~ mean misfit angle, Eq. 37).
#
# Returns: list(y, mu, wt_final, eigvals, eigenvalue_gap, n_iter)
# ===
.wissi_stage2 <- function(normals, slips, weights, y_init,
                          max_iter = 50L, tol_deg = 1e-4) {
  N <- nrow(normals)
  lambda_max <- sqrt(3) / 2
  y_cur <- y_init
  mu_i <- rep(1, N)

  for (iter in seq_len(max_iter)) {
    sigma <- .y6_to_sigma(y_cur)

    # Per-fault magnitude weight: tau_mag / lambda_max in (0, 1]
    mu_i <- vapply(seq_len(N), function(i) {
      max(1e-4, .tau_mag(sigma, normals[i, ]) / lambda_max)
    }, numeric(1L))

    wt_new <- weights * mu_i
    M5_new <- .build_M5(normals, slips, wt_new)
    res <- .solve_M5(M5_new)
    y_new <- .resolve_sense(res$y, normals, slips)

    # Convergence: Theta between successive iterates
    delta <- angular_stress_distance(y_cur, y_new)
    y_cur <- y_new
    if (delta < tol_deg) break
  }

  gap <- res$eigvals[4L] - res$eigvals[5L]
  list(
    y = y_cur, mu = mu_i, wt_final = wt_new,
    eigvals = res$eigvals, eigenvalue_gap = gap, n_iter = iter
  )
}


# ===
# STAGE 3 — Sense annealing (robust slip sense handling)
# ===
# Replaces the binary majority-vote sense resolution with a continuous
# tanh-weighted annealing schedule. At low gamma (early iterations) the method
# is sense-agnostic (robust, like Hansen 2013). As gamma increases, it
# commits progressively to the predicted sense (accurate, like Angelier 1990).
#
# For each fault:
#   phi_i = tanh(gamma * tau_i . slip_i)  in (-1, 1)
#   - phi_i > 0: slip sense consistent with current tensor  -> normal weight
#   - phi_i < 0: slip sense inconsistent                    -> flip slip
#   - phi_i ~ 0: ambiguous                                  -> low weight
#
# Combined weight: omega_i * mu_i * |phi_i|
# Slips are flipped where phi_i < 0 before building M5.
#
# Returns: list(y, phi, mu, slips_corrected, n_flipped, wt_final,
#               eigvals, eigenvalue_gap, n_iter_total)
# ===
.wissi_stage3 <- function(normals, slips, weights, y_init,
                          gamma_max = 10,
                          n_anneal = 8L,
                          max_iter = 30L,
                          tol_deg = 1e-4) {
  N <- nrow(normals)
  lambda_max <- sqrt(3) / 2
  y_cur <- y_init
  phi_i <- rep(1, N)
  mu_i <- rep(1, N)
  slips_adj <- slips
  n_iter_tot <- 0L

  for (ann in seq_len(n_anneal)) {
    gamma <- gamma_max * (ann / n_anneal) # linear schedule 0 -> gamma_max

    for (iter in seq_len(max_iter)) {
      sigma <- .y6_to_sigma(y_cur)

      # Soft sense weights
      tau_dot_s <- vapply(seq_len(N), function(i) {
        sum(.shear_tau(sigma, normals[i, ]) * slips[i, ])
      }, numeric(1L))
      phi_i <- tanh(gamma * tau_dot_s)

      # Magnitude weights
      mu_i <- vapply(seq_len(N), function(i) {
        max(1e-4, .tau_mag(sigma, normals[i, ]) / lambda_max)
      }, numeric(1L))

      # Adjust slip sense: flip where phi_i < 0
      slips_adj <- slips
      slips_adj[phi_i < 0, ] <- -slips[phi_i < 0, ]

      # Combined weight: quality * magnitude * |sense confidence|
      wt_i <- weights * mu_i * pmax(abs(phi_i), 0.01)

      M5_new <- .build_M5(normals, slips_adj, wt_i)
      res <- .solve_M5(M5_new)
      y_new <- .resolve_sense(res$y, normals, slips_adj)

      delta <- angular_stress_distance(y_cur, y_new)
      n_iter_tot <- n_iter_tot + 1L
      y_cur <- y_new
      if (delta < tol_deg) break
    }
  }

  gap <- res$eigvals[4L] - res$eigvals[5L]
  list(
    y = y_cur,
    phi = phi_i,
    mu = mu_i,
    slips_corrected = slips_adj,
    n_flipped = sum(phi_i < 0),
    wt_final = wt_i,
    eigvals = res$eigvals,
    eigenvalue_gap = gap,
    n_iter_total = n_iter_tot
  )
}


# ===
# STAGE 4 — Analytic uncertainty via eigenvalue perturbation theory
# ===
# Propagates slip direction measurement noise (sigma_alpha_deg) through the
# eigenproblem to obtain Cov(x5) in 5D sigma-space, then maps to Cov(y6).
#
# For each fault, the derivative of M5 w.r.t. a unit angular perturbation of
# the slip direction is computed by finite difference (delta = 1e-5 rad).
# The perturbation covariance is:
#
#   Cov(x5) ≈ sigma_alpha^2 * sum_i sum_{j!=min} [v_j^T dM5_i x5]^2
#                                                 / (lam_min - lam_j)^2
#                                                 * v_j v_j^T
#
# The eigenvalue gap (lam_2 - lam_1) serves as a condition number:
# a small gap means the solution direction is poorly separated from the
# next eigenvector, and uncertainty estimates become unreliable.
#
# Returns: list(Cov5, Cov_y6, eigval_gap, cov_eigvals,
#               sigma1_unc_deg, Phi_unc)
# ===
.wissi_stage4 <- function(normals, slips, weights, y_opt,
                          sigma_alpha_deg = 10) {
  N <- nrow(normals)
  sigma_alpha <- sigma_alpha_deg * pi / 180

  M5 <- .build_M5(normals, slips, weights)
  eig <- eigen(M5, symmetric = TRUE)
  lam <- eig$values # decreasing: lam[1] >= ... >= lam[5]
  V <- eig$vectors

  gap <- lam[4L] - lam[5L] # second - smallest eigenvalue gap

  x5 <- V[, 5L] # optimal in 5D
  Cov5 <- matrix(0, 5, 5)

  for (i in seq_len(N)) {
    n <- normals[i, ]
    s <- slips[i, ]
    ep6 <- .eps6_prime(s, n)

    # In-plane perturbation direction: b = n x s
    b <- c(n[2] * s[3] - n[3] * s[2], n[3] * s[1] - n[1] * s[3], n[1] * s[2] - n[2] * s[1])
    b_norm <- sqrt(sum(b^2))
    if (b_norm < 1e-12) next
    b <- b / b_norm

    # Numerical derivative of ep5 w.r.t. slip rotation toward b
    delta_a <- 1e-5
    s_pert <- cos(delta_a) * s + sin(delta_a) * b
    s_pert <- s_pert - sum(s_pert * n) * n
    s_pert <- s_pert / sqrt(sum(s_pert^2))

    dep6 <- (.eps6_prime(s_pert, n) - ep6) / delta_a
    dep5 <- drop(crossprod(.B, dep6))

    # dM5 = weights[i] * (dep5 ep5^T + ep5 dep5^T)
    ep5 <- drop(crossprod(.B, ep6))
    dM5 <- weights[i] * (tcrossprod(dep5, ep5) + tcrossprod(ep5, dep5))

    # Perturbation theory: contributions from non-minimum eigenvectors
    for (j in 1:4) {
      vj <- V[, j]
      num <- drop(t(vj) %*% dM5 %*% x5)
      denom <- lam[5L] - lam[j] # negative
      if (abs(denom) > 1e-12) {
        Cov5 <- Cov5 + sigma_alpha^2 * (num / denom)^2 * tcrossprod(vj)
      }
    }
  }

  Cov_y6 <- .B %*% Cov5 %*% t(.B)
  cov_eig <- eigen(Cov5, symmetric = TRUE)$values
  # 1-sigma uncertainty on principal stress axis orientations (approximate)
  sigma1_unc <- sqrt(pmax(0, max(cov_eig))) * 180 / pi
  # Phi uncertainty: propagate through eigenvalue difference ratio
  # approximate as the RMS of diagonal Cov5 mapped to Phi sensitivity
  Phi_unc <- sqrt(sum(pmax(0, diag(Cov5)))) * 180 / pi / 90

  list(
    Cov5 = Cov5,
    Cov_y6 = Cov_y6,
    eigval_gap = gap,
    cov_eigvals = cov_eig,
    sigma1_unc_deg = sigma1_unc,
    Phi_unc = Phi_unc
  )
}


# ===
# STAGE 5 — Polyphase separation via spectral clustering on S^5
# ===
# Each fault is represented by the normalised pole of its solution great
# hypercircle in 5D (the ep5 vector). Faults from the same stress phase
# tend to have similar poles. The algorithm:
#
#   1. Compute pairwise angular distances between fault poles.
#   2. Build the Gaussian affinity matrix \eqn{K_ij = \exp(-D_{ij}^2 / 2*\sigma_K^2)}.
#   3. Compute the normalised graph Laplacian L_sym.
#   4. Select k via the eigengap heuristic on the k_max smallest eigenvalues.
#   5. k-means cluster the row-normalised spectral embedding.
#   6. Run wissi() on each cluster to obtain per-phase tensors.
#
# sigma_K_deg controls the affinity bandwidth — faults within sigma_K_deg
# of each other (in the ASD sense) are considered similar. A value of 20-40
# degrees is typically appropriate; larger values merge nearby clusters.
#
# Returns: list(assignment, k_opt, gaps, phase_results, D_mat)
#   assignment    : integer vector of length N (phase label per fault)
#   k_opt         : number of phases identified
#   gaps          : Laplacian eigenvalue gaps (eigengap criterion)
#   phase_results : list of k wissi() results, one per phase
#   D_mat         : N x N pairwise ASD matrix between fault poles
# ===
#' @importFrom stats kmeans
.wissi_stage5 <- function(normals, slips, weights,
                          k_max = 4L,
                          sigma_K_deg = 30,
                          wissi_args = list(),
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(normals)

  # Fault pole vectors in 5D (normalised ep5)
  poles <- matrix(0, N, 5)
  for (i in seq_len(N)) {
    ep5 <- drop(crossprod(.B, .eps6_prime(slips[i, ], normals[i, ])))
    poles[i, ] <- ep5 / sqrt(sum(ep5^2))
  }

  # Pairwise ASD between poles (abs for antipodal equivalence)
  D_mat <- matrix(0, N, N)
  for (i in seq_len(N - 1L)) {
    for (j in (i + 1L):N) {
      d <- acos(pmax(-1, pmin(1, abs(sum(poles[i, ] * poles[j, ]))))) * 180 / pi
      D_mat[i, j] <- D_mat[j, i] <- d
    }
  }

  # Gaussian affinity matrix
  K <- exp(-D_mat^2 / (2 * sigma_K_deg^2))

  # Normalised symmetric Laplacian: L = I - D^{-1/2} K D^{-1/2}
  d_inv_sqrt <- 1 / sqrt(pmax(rowSums(K), 1e-12))
  L_sym <- diag(N) - d_inv_sqrt * K * rep(d_inv_sqrt, each = N)

  # Eigenvectors of L_sym (smallest eigenvalues carry cluster structure)
  eig_L <- eigen(L_sym, symmetric = TRUE)
  lam_inc <- rev(eig_L$values) # increasing order
  gaps <- diff(lam_inc[seq_len(k_max + 1L)])
  k_opt <- max(1L, which.max(gaps))

  # Spectral embedding: k_opt smallest eigenvectors (from end of decreasing)
  U <- eig_L$vectors[, (N - k_opt + 1L):N, drop = FALSE]
  r_norm <- sqrt(rowSums(U^2))
  r_norm[r_norm < 1e-10] <- 1
  U_norm <- U / r_norm

  # k-means clustering in the spectral embedding
  km <- stats::kmeans(U_norm, centers = k_opt, nstart = 25L, iter.max = 200L)
  assignment <- km$cluster

  # Run wissi() on each cluster
  phase_results <- vector("list", k_opt)
  for (ph in seq_len(k_opt)) {
    idx <- which(assignment == ph)
    if (length(idx) < 4L) {
      phase_results[[ph]] <- list(error = "fewer than 4 faults in phase")
      next
    }
    args <- c(
      list(
        normals = normals[idx, , drop = FALSE],
        slips = slips[idx, , drop = FALSE],
        weights = weights[idx]
      ),
      wissi_args
    )
    phase_results[[ph]] <- tryCatch(do.call(wissi, args),
      error = function(e) list(error = e$message)
    )
  }

  list(
    assignment = assignment,
    k_opt = k_opt,
    gaps = gaps,
    phase_results = phase_results,
    D_mat = D_mat
  )
}


# ===
# WISSI — Main user-facing function
# ===
#
#' Weighted Iterative Sigma-Space Inversion (WISSI)
#'
#' @description
#' Combines four classic fault slip inversion algorithms into a single
#' coherent framework operating in the Yamaji-Sato 5-sphere sigma-space.
#' Integrates the linear eigenvector approach of Yamaji & Sato (2006),
#' the magnitude-aware upsilon criterion of Angelier (1990) and Mostafa (2005),
#' the sense annealing strategy inspired by Hansen (2013), and the analytic
#' uncertainty framework of Michael (1984), extended to the curved geometry
#' of the 5-sphere.
#'
#' @inheritParams slip_inversion
#' @param normals \code{N x 3} numeric matrix of unit fault plane normals in
#'   a right-handed Cartesian coordinate system (x = East, y = North, z = Up).
#'   Normals should point toward the footwall (upward-pointing for
#'   sub-horizontal faults).
#' @param slips \code{N x 3} numeric matrix of unit slip direction vectors,
#'   parallel to the shear traction and pointing in the direction of
#'   hanging-wall motion. Must be perpendicular to the corresponding row of
#'   \code{normals}.
#' @param weights Optional length-\code{N} numeric vector of non-negative data
#'   quality weights. \code{NA} values are replaced with the observed mean
#'   before scaling (neutral imputation). The vector is normalised internally
#'   so that \code{mean(weights) = 1}, meaning only the ratios between weights
#'   matter. Use [scale_weights()] to construct weights from field quality
#'   ranks, measurement errors, and prior RUP values. Default: uniform
#'   weights.
#' @param sigma_alpha_deg Estimated angular measurement error of slip
#'   directions in degrees. Used in Stage 4 (analytic uncertainty) to scale
#'   the perturbation covariance matrix. A value of 10 &deg; is a
#'   conservative default for well-measured slickenlines; increase to 20-30
#'   &deg; for less certain measurements. Default: \code{10}.
#' @param gamma_max Maximum sharpness parameter for the sense annealing
#'   schedule (Stage 3). Controls how strongly the final solution commits to
#'   the predicted slip sense. Higher values produce sharper sense
#'   discrimination; set to \code{0} to disable sense annealing entirely,
#'   making the inversion fully sense-agnostic in the style of Hansen (2013).
#'   Default: \code{10}.
#' @param n_anneal Number of annealing steps in the outer sense annealing loop
#'   (Stage 3). The sharpness parameter \code{gamma} increases linearly from
#'   \code{0} to \code{gamma_max} over \code{n_anneal} steps. More steps give
#'   a smoother transition between sense-agnostic and sense-committed
#'   behaviour. Default: \code{8}.
#' @param max_iter Maximum number of inner iterations per annealing step.
#'   Default: \code{50}.
#' @param tol_deg Convergence tolerance in angular stress distance (degrees).
#'   The inner loop terminates when the angular stress distance between
#'   successive iterates falls below this value. Because angular stress
#'   distance approximates the mean misfit angle (Yamaji & Sato 2006, Eq. 37),
#'   this tolerance has a direct physical interpretation. Default: \code{1e-4}.
#' @param run_stage4 Logical. If \code{TRUE}, compute the analytic uncertainty
#'   estimates via eigenvalue perturbation theory (Stage 4) and return them in
#'   the \code{unc} element of the result. Setting to \code{FALSE} skips this
#'   step and is useful for bootstrap resampling where many inversions are
#'   performed. Default: \code{TRUE}.
#' @param robust Logical. If \code{TRUE}, apply iteratively reweighted robust
#'   estimation using the kernel specified by \code{robust_kernel}. Faults
#'   with large angular misfits are progressively downweighted at each
#'   iteration, reducing the influence of outliers on the final tensor.
#'   Default: \code{TRUE}.
#' @param robust_kernel Character string specifying the robust weighting
#'   kernel. One of:
#'   \describe{
#'     \item{\code{"huber"}}{Huber kernel: unit weight below the threshold
#'       \code{robust_k_deg}, declining as \code{k/alpha} above it. Gradual
#'       downweighting; the recommended default for most datasets.}
#'     \item{\code{"tukey"}}{Tukey bisquare kernel: unit weight near zero
#'       misfit, declining smoothly to exactly zero at \code{robust_k_deg}.
#'       Faults beyond the threshold are completely excluded. Use when
#'       genuinely bad measurements are known to be present.}
#'     \item{\code{"andrews"}}{Andrews sine kernel: similar to Tukey but
#'       differentiable at the threshold. Slightly smoother rejection profile.}
#'   }
#'   Default: \code{"huber"}.
#' @param robust_k_deg Threshold angle in degrees for the robust kernel.
#'   Faults with angular misfit below this value receive full (or near-full)
#'   weight; faults above it are progressively downweighted or excluded
#'   depending on the chosen kernel. Should be set relative to the expected
#'   data noise level: \code{30}&deg; is appropriate for typical
#'   slickenline data, \code{20} &deg; for high-quality datasets, and \code{45} &deg; when
#'   significant scatter is expected. Default: \code{30}.
#' @param robust_mad_k Multiplier for the MAD-based outlier flagging applied
#'   after convergence. Faults whose angular misfit exceeds
#'   \code{median(alpha) + robust_mad_k * 1.4826 * MAD(alpha)} are reported
#'   in \code{outliers_mad}. This is a diagnostic flag only and does not
#'   affect the inversion result. Smaller values flag more faults; \code{2.5}
#'   corresponds roughly to the 99th percentile under a normal distribution.
#'   Default: \code{2.5}.
#' @param init_median Logical. If \code{TRUE} and \code{robust = TRUE},
#'   initialise the iteration from the spatial median on S\eqn{^5} computed
#'   via the Weiszfeld algorithm, rather than the standard weighted eigenvector.
#'   The spatial median minimises the sum of angular stress distances to all
#'   single-fault solution poles and is substantially more resistant to
#'   outliers than the mean-based eigenvector initialisation. Has no effect
#'   when \code{robust = FALSE}. Default: \code{TRUE}.
#'
#' @details
#' ## Algorithm stages
#'
#' *WISSI* proceeds through five sequential stages, each addressing a specific
#' limitation of the classical inversion methods.
#'
#' **Stage 1 — Initialisation.** Constructs the weighted moment matrix
#' \eqn{M_5 = \sum_i \omega_i \epsilon'_i \epsilon'^T_i} in the 5D deviator
#' subspace of the Yamaji-Sato sigma-space and returns its
#' smallest-eigenvalue eigenvector as the initial stress estimate. When
#' \code{robust = TRUE} and \code{init_median = TRUE}, the Weiszfeld
#' algorithm on \eqn{S^5} is used instead, providing a starting point
#' that is resistant to leverage points from badly oriented or erroneous
#' fault data.
#'
#' **Stage 2 — Magnitude reweighting.** Translates the Angelier (1990) and
#' Mostafa (2005) upsilon (\eqn{\upsilon_i}) magnitude criterion into the
#' sigma-space eigenproblem by replacing the global \eqn{\lambda = \sqrt{3}/2}
#' with per-fault weights \eqn{\mu_i = |\tau_i| / \lambda_{\max}}.
#' Mechanically degenerate fault planes (those containing a principal stress
#' axis, which produce near-zero shear traction) are automatically
#' downweighted without requiring any explicit degeneracy test.
#' Convergence is measured by the angular stress distance \eqn{\Theta} between
#' successive iterates, which approximates the mean misfit angle
#' (Yamaji & Sato 2006, Eq. 37) and therefore has a direct physical
#' interpretation.
#'
#' **Stage 3 — Sense annealing.** Replaces the binary majority-vote sense
#' resolution of Yamaji & Sato (2006) with a continuous
#' \eqn{\tanh(\gamma \cdot \hat{\tau}_i \cdot \hat{s}_i)} weighting schedule.
#' At small \eqn{\gamma} (early annealing steps) the method is effectively
#' sense-agnostic, as in Hansen (2013), providing robustness when slip senses
#' are uncertain or partially incorrect. As \eqn{\gamma} increases toward
#' \code{gamma_max}, the method commits progressively to the predicted slip
#' sense, recovering the accuracy of Angelier (1990) for reliable data.
#' Faults whose predicted and recorded senses are inconsistent have their
#' slip vectors flipped internally; the number of such corrections is reported
#' in \code{n_flipped_sense}.
#'
#' **Stage 4 — Analytic uncertainty.** Propagates slip direction measurement
#' noise (\code{sigma_alpha_deg}) through the eigenproblem via first-order
#' perturbation theory, yielding a \eqn{5 \times 5} covariance matrix in
#' sigma-space and approximate \eqn{1\sigma} confidence bounds on the
#' principal stress axis orientations and \eqn{\Phi}. This extends the
#' analytic covariance framework of Michael (1984) to the geometrically
#' correct curved space of the 5-sphere. The eigenvalue gap
#' \eqn{\lambda_2 - \lambda_1} of \eqn{M_5} serves as a condition number:
#' a small gap indicates that the solution is poorly separated from the next
#' eigenvector and that uncertainty estimates may be unreliable.
#'
#' **Stage 5 (polyphase) — Spectral clustering.** Available via
#' \code{\link{slip_inversion_wissi_polyphase}}. Represents each fault by the pole of its
#' solution great hypercircle on \eqn{S^5}, builds a Gaussian affinity matrix
#' using the angular stress distance \eqn{\Theta}, and applies normalised
#' spectral clustering with automatic phase count selection via the eigengap
#' heuristic. WISSI is then run independently on each identified phase.
#'
#' ## Robust estimation
#'
#' When \code{robust = TRUE}, an additional misfit-based weight
#' \eqn{w_i^{\text{rob}} = \kappa(\alpha_i)} is computed at each iteration
#' and multiplied into the combined weight \eqn{\tilde{\omega}_i = \omega_i
#' \cdot \mu_i \cdot |\phi_i| \cdot w_i^{\text{rob}}}, where \eqn{\alpha_i}
#' is the unsigned angular misfit between the predicted shear traction and the
#' observed slip direction, \eqn{\mu_i} is the magnitude weight from Stage 2,
#' and \eqn{\phi_i} is the sense confidence from Stage 3. The three available
#' kernels \eqn{\kappa} differ in tail behaviour: the Huber kernel applies a
#' gradual \eqn{k/\alpha} decay beyond the threshold and is the recommended
#' default; the Tukey bisquare kernel hard-zeros faults beyond the threshold
#' for aggressive rejection of confirmed bad measurements; the Andrews sine
#' kernel provides a smooth differentiable alternative to Tukey. After
#' convergence, faults are additionally screened using a
#' median-absolute-deviation (MAD) criterion that adapts automatically to the
#' noise level of the specific dataset, reported in \code{outliers_mad}. This
#' flag is purely diagnostic and does not alter the inversion result.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{sigma}}{3x3 reduced stress tensor in the input Cartesian
#'     frame.}
#'   \item{\code{y}}{6D unit y-vector on \eqn{S^5} representing the optimal
#'     stress tensor.}
#'   \item{\code{sigma1}, \code{sigma2}, \code{sigma3}}{Unit vectors of the
#'     principal stress axes, from maximum (\eqn{\sigma_1}) to minimum
#'     (\eqn{\sigma_3}) compressive stress.}
#'   \item{\code{eigenvalues}}{Eigenvalues of \code{sigma} in decreasing
#'     order.}
#'   \item{`principal_axes`}{unit vectors of principal stress axes (max to min)}
#'   \item{`principal_vals`}{eigenvalues of `stress_tensor` (decreasing)}
#'   \item{\code{Phi}}{Shape ratio
#'     \eqn{\Phi = (\sigma_1 - \sigma_2)/(\sigma_1 - \sigma_3) \in [0, 1]}.}
#'   \item{\code{alpha_deg}}{Per-fault unsigned angular misfit in degrees
#'     (0-90), the standard line metric.}
#'   \item{\code{alpha_signed_deg}}{Per-fault signed angular misfit in degrees
#'     (0-180). Values above 90 indicate that the recorded slip sense is
#'     opposite to the predicted shear traction direction.}
#'   \item{\code{mean_alpha}}{Mean unsigned angular misfit across all faults
#'     (degrees).}
#'   \item{\code{median_alpha}}{Median unsigned angular misfit (degrees). More
#'     robust than the mean in the presence of outliers.}
#'   \item{\code{suspected_flipped}}{Integer vector of row indices where
#'     \code{alpha_signed_deg > 90}, flagging faults whose recorded slip sense
#'     may be incorrect.}
#'   \item{\code{outliers_mad}}{Integer vector of row indices flagged as
#'     outliers by the MAD criterion (only when \code{robust = TRUE}).}
#'   \item{\code{n_flipped_sense}}{Number of faults whose slip sense was
#'     corrected internally during Stage 3.}
#'   \item{\code{slips_corrected}}{N x 3 matrix of sense-corrected slip
#'     vectors used in the final inversion.}
#'   \item{\code{mu}}{Per-fault magnitude weights \eqn{\mu_i} from Stage 2.}
#'   \item{\code{phi_sense}}{Per-fault \eqn{\tanh} sense confidence values
#'     from Stage 3. Positive values indicate consistent sense; negative
#'     values indicate corrected (flipped) senses.}
#'   \item{\code{w_robust}}{Per-fault robust kernel weights \eqn{w_i^{\text{rob}}}
#'     from the final iteration (only when \code{robust = TRUE}).}
#'   \item{\code{wt_combined}}{Per-fault combined weights
#'     \eqn{\tilde{\omega}_i} used in the final M5 matrix.}
#'   \item{\code{eigenvalue_gap}}{Difference \eqn{\lambda_2 - \lambda_1} of
#'     the final \eqn{M_5} matrix. Larger values indicate a better-conditioned
#'     solution.}
#'   \item{\code{M5_eigvals}}{All five eigenvalues of the final \eqn{M_5}
#'     matrix in decreasing order.}
#'   \item{\code{unc}}{List of analytic uncertainty estimates from Stage 4
#'     (only when \code{run_stage4 = TRUE}), containing \code{Cov5},
#'     \code{Cov_y6}, \code{eigval_gap}, \code{cov_eigvals},
#'     \code{sigma1_unc_deg}, and \code{Phi_unc}.}
#'   \item{\code{n_iter_total}}{Total number of inner iterations performed
#'     across all annealing steps.}
#' }
#'
#'
#' @references
#'   Angelier, J. (1990). Inversion of field data in fault tectonics to obtain
#'   the regional stress. III. A new rapid direct inversion method by
#'   analytical means. \emph{Geophys. J. Int.}, 103, 363-376.
#'
#'   Hansen, J.A. (2013). Direct inversion of stress from sets of fault slip
#'   data with unknown slip sense. \emph{J. Struct. Geol.}, 51, 54-65.
#'
#'   Michael, A.J. (1984). Determination of stress from slip data: faults and
#'   folds. \emph{J. Geophys. Res.}, 89, 11517-11526.
#'
#'   Mostafa, M.E. (2005). Iterative linear inversion of stress tensor from
#'   fault-slip data. \emph{J. Struct. Geol.}, 27, 930-940.
#'
#'   Pascal, C. (2022). \emph{Paleostress Inversion Techniques}. Elsevier.
#'
#'   Yamaji, A. & Sato, K. (2006). Distances for the solutions of stress
#'   tensor inversion in relation to misfit angles that accompany the
#'   solutions. \emph{Geophys. J. Int.}, 167, 933-942.
#'   
#'   Stephan (in prep)
#'   
#' @name slip_inversion_wissi
#' 
#' @family stress-inversion
#' @family wissi
#'
#' @examples
#' set.seed(20250411)
#' 
#' nx <- length(angelier1990)
#' par(mfrow = c(2, nx/2))
#'
#' invisible(lapply(seq_len(nx), function(i) {
#'   # inversion
#'   x <- angelier1990[[i]]
#'   res <- slip_inversion_wissi(x)
#'
#'   # some stress shape
#'   phi_val <- round(res$stress_shape$phi, 2)
#'
#'   # misfit
#'   rup_val <- round(res$misfit$rup_mean, 2)
#'
#'   # Plot the faults (color-coded by RUP%) and show the principal stress axes
#'   stereoplot(title = names(angelier1990)[i], guides = FALSE)
#'   stereo_shmax(res$SHmax)
#'   fault_plot(x, col = assign_col(res$misfit$rup))
#'   points(res$principal_axes, col = 1:3, pch = 16, cex = 1.5)
#'   text(res$principal_axes,
#'     label = rownames(res$principal_axes),
#'     col = 1:3, adj = -.25
#'   )
#'   legend("topleft", col = 2:4, legend = rownames(res$principal_axes), pch = 16)
#'   title(sub = bquote(Phi == .(phi_val) ~ "|" ~ bar("RUP") == .(rup_val) * "%"))
#' }))
NULL

#' @rdname slip_inversion_wissi
#' @noRd
wissi <- function(normals,
                  slips,
                  weights = NULL,
                  sigma_alpha_deg = 10,
                  gamma_max = 10,
                  n_anneal = 8L,
                  max_iter = 50L,
                  tol_deg = 1e-4,
                  run_stage4 = TRUE,
                  # --- new robust options ---
                  robust          = TRUE,
                  robust_kernel   = c("huber", "tukey", "andrews"),
                  robust_k_deg    = 30,
                  robust_mad_k    = 2.5,
                  init_median     = TRUE
                  ) {
  # --- Input validation ---
  if (!is.matrix(normals) || !is.matrix(slips)) {
    stop("'normals' and 'slips' must be matrices.")
  }
  if (ncol(normals) != 3L || ncol(slips) != 3L) {
    stop("'normals' and 'slips' must have 3 columns (x, y, z).")
  }
  if (nrow(normals) != nrow(slips)) {
    stop("'normals' and 'slips' must have the same number of rows.")
  }
  if (nrow(normals) < 4L) {
    stop("At least 4 fault slip measurements are required.")
  }

  robust_kernel <- match.arg(robust_kernel)
  # Kernel weight function
  .w_robust <- switch(robust_kernel,
                      huber    = function(a) .w_huber(a,   k = robust_k_deg),
                      tukey    = function(a) .w_tukey(a,   k = robust_k_deg),
                      andrews  = function(a) .w_andrews(a, k = robust_k_deg)
  )
  
  N <- nrow(normals)

  # --- Weights ---
  if (is.null(weights)) {
    wt <- rep(1, N)
  } else {
    stopifnot(
      is.numeric(weights), length(weights) == N,
      any(!is.na(weights)), all(weights >= 0, na.rm = TRUE)
    )
    weights[is.na(weights)] <- mean(weights, na.rm = TRUE)
    wt <- weights / mean(weights)
  }

  # --- Stage 1: initialisation ---
  # Robust: spatial median on S5; otherwise: standard weighted eigenvector
  if (robust && init_median) {
    x5_init <- .s5_median(normals, slips, wt)
    y_init  <- drop(.B %*% x5_init)
    y_init  <- .resolve_sense(y_init / sqrt(sum(y_init^2)), normals, slips)
  } else {
    s1     <- .wissi_stage1(normals, slips, wt)
    y_init <- s1$y
  }
  
  # --- Stages 2 + 3 combined: magnitude + sense + robust reweighting ---
  lambda_max <- sqrt(3) / 2
  y_cur      <- y_init
  phi_i      <- rep(1, N)
  mu_i       <- rep(1, N)
  w_rob      <- rep(1, N)
  slips_adj  <- slips
  n_iter_tot <- 0L
  
  for (ann in seq_len(n_anneal)) {
    gamma <- gamma_max * (ann / n_anneal)
    
    for (iter in seq_len(max_iter)) {
      sigma <- .y6_to_sigma(y_cur)
      
      # --- Magnitude weights (Stage 2 / Mostafa) ---
      mu_i <- vapply(seq_len(N), function(i)
        max(1e-4, .tau_mag(sigma, normals[i,]) / lambda_max), numeric(1L))
      
      # --- Sense annealing (Stage 3) ---
      tau_dot_s <- vapply(seq_len(N), function(i)
        sum(.shear_tau(sigma, normals[i,]) * slips[i,]), numeric(1L))
      phi_i     <- tanh(gamma * tau_dot_s)
      slips_adj <- slips
      slips_adj[phi_i < 0, ] <- -slips[phi_i < 0, ]
      
      # --- Robust misfit weights (new) ---
      if (robust) {
        alpha_i <- vapply(seq_len(N), function(i) {
          tau <- .shear_tau(sigma, normals[i,])
          tn  <- sqrt(sum(tau^2))
          if (tn < 1e-12) return(90)
          acos(pmax(-1, pmin(1, abs(sum(tau/tn * slips_adj[i,]))))) * 180/pi
        }, numeric(1L))
        w_rob <- .w_robust(alpha_i)
      }
      
      # Combined weight: quality * magnitude * |sense| * robust
      wt_i <- wt * mu_i * pmax(abs(phi_i), 0.01) * w_rob
      
      M5_new <- .build_M5(normals, slips_adj, wt_i)
      res    <- .solve_M5(M5_new)
      y_new  <- .resolve_sense(res$y, normals, slips_adj)
      
      delta      <- angular_stress_distance(y_cur, y_new)
      n_iter_tot <- n_iter_tot + 1L
      y_cur      <- y_new
      if (delta < tol_deg) break
    }
  }
  
  # --- Stage 4: analytic uncertainty ---
  unc <- if (run_stage4){
    .wissi_stage4(normals, slips_adj, wt_i, y_cur,
                  sigma_alpha_deg = sigma_alpha_deg)
  } else NULL
  
  # --- Post-processing ---
  dec      <- .decompose_y(y_cur)
  sigma    <- dec$sigma
  a_deg    <- .alpha_deg(sigma, normals, slips_adj)
  a_signed <- .alpha_signed_deg(sigma, normals, slips_adj)
  
  # MAD-based outlier flags
  outlier_mad <- if (robust){
    .flag_outliers_mad(a_deg, k = robust_mad_k)
  }  else rep(FALSE, N)
  
  list(
    sigma             = sigma,
    y                 = y_cur,
    sigma1            = dec$sigma1,
    sigma2            = dec$sigma2,
    sigma3            = dec$sigma3,
    eigenvalues       = dec$eigenvalues,
    Phi               = dec$Phi,
    alpha_deg         = a_deg,
    alpha_signed_deg  = a_signed,
    mean_alpha        = mean(a_deg, na.rm = TRUE),
    median_alpha      = median(a_deg, na.rm = TRUE),
    suspected_flipped = which(a_signed > 90),
    outliers_mad      = which(outlier_mad),       # MAD-flagged outliers
    n_flipped_sense   = sum(phi_i < 0),
    slips_corrected   = slips_adj,
    mu                = mu_i,
    phi_sense         = phi_i,
    w_robust          = w_rob,                    # final robust weights
    wt_combined       = wt_i,                     # all weights combined
    eigenvalue_gap    = res$eigvals[4L] - res$eigvals[5L],
    M5_eigvals        = res$eigvals,
    unc               = unc,
    n_iter_total      = n_iter_tot
  )
}

#' @rdname slip_inversion_wissi
#' @export
slip_inversion_wissi <- function(x, weights = NULL,
                                 sigma_alpha_deg = 10,
                                 gamma_max = 10,
                                 n_anneal = 8L,
                                 max_iter = 50L,
                                 tol_deg = 1e-4,
                                 run_stage4 = TRUE,
                                 # --- new robust options ---
                                 robust          = TRUE,
                                 robust_kernel   = c("huber", "tukey", "andrews"),
                                 robust_k_deg    = 30,
                                 robust_mad_k    = 2.5,
                                 init_median     = TRUE) {
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  res <- wissi(normals, slips, weights = weights, sigma_alpha_deg = sigma_alpha_deg,
               gamma_max = gamma_max, n_anneal = n_anneal, max_iter = max_iter,
               tol_deg = tol_deg, run_stage4 = run_stage4, 
               robust = robust, robust_kernel = robust_kernel, robust_k_deg = robust_k_deg,
               robust_mad_k = robust_mad_k, init_median = init_median)
  
  #res$principal_axes <- as.Vec3(rbind(res$sigma1, res$sigma2, res$sigma3)) |> Line()
  #res$principal_vals <- res$eigenvalues
  
  s <- sigma2stress(res$sigma)
  res$principal_axes <- s$axes
  res$principal_vals <- s$vals
  res$sigma1 <- res$sigma2 <- res$sigma3 <- res$eigenvalues <- res$Phi <- NULL
  
  res$stress_shape <- stress_shape(res$sigma)
  
  res$SHmax <- SH(res$principal_axes[1, ], res$principal_axes[2, ], res$principal_axes[3, ], R = res$stress_shape$R)
  
  res$slips_corrected <- as.Vec3(res$slips_corrected)
  
  res$misfit <- slip_inversion_misfit(res$sigma, x)
  
  return(res)
}


# WISSI_POLYPHASE — Polyphase separation wrapper
wissi_polyphase <- function(normals, slips,
                            weights = NULL,
                            k_max = 4L,
                            sigma_K_deg = 30,
                            seed = NULL,
                            ...) {
  
  if (!is.matrix(normals) || !is.matrix(slips)) {
   stop("'normals' and 'slips' must be matrices.")
  }
  if (nrow(normals) < 4L) {
    stop("At least 4 fault slip measurements are required.")
  }
  
  N <- nrow(normals)
  wt <- if (is.null(weights)) {
    rep(1, N)
  } else {
    w <- weights
    w[is.na(w)] <- mean(w, na.rm = TRUE)
    w / mean(w)
  }
  
  .wissi_stage5(normals, slips, wt,
                k_max       = k_max,
                sigma_K_deg = sigma_K_deg,
                wissi_args  = list(...),
                seed        = seed
  )
}

#' Polyphase stress inversion via spectral clustering on \eqn{S^5} (Stage 5).
#'
#' Identifies `k` stress phases automatically using the *eigengap heuristic*,
#' then runs `slip_inversion_wissi()` on each phase subset.
#'
#' @inheritParams slip_inversion_wissi
#' @param k_max       Maximum number of phases to consider. Default `4`.
#' @param sigma_K_deg Affinity bandwidth in degrees of angular stress
#'                    distance. Faults within this distance are considered
#'                    similar. Default `30.` Increase to merge nearby phases,
#'                    decrease to split them.
#' @param seed        Optional RNG seed for k-means reproducibility.
#' @param ...         Additional arguments passed to `slip_inversion_wissi()` for each phase.
#'
#' @return A named list with: \describe{
#'   \item{`assignment`}{integer vector of length N (phase label 1..k per fault)}
#'   \item{`k_opt`}{number of phases identified}
#'   \item{`gaps`}{Laplacian eigenvalue gaps (*eigengap* criterion)}
#'   \item{`phase_results`}{list of `k` `slip_inversion_wissi()` results, one per phase}
#'   \item{`D_mat`}{N x N pairwise ASD matrix between fault poles (degrees)}
#'   }
#'  
#' @family wissi 
#' 
#' @export
#' 
#' @examples
#' res <- slip_inversion_wissi_polyphase(angelier1990$KAM, sigma_K_deg = 15)
#' 
#' # check amount of clusters detected
#' res$k_opt
#' 
#' 
#' # plot the results in lambert projection:
#' cols <- assign_col_d(seq_len(res$k_opt))#' 
#' par(mfrow = c(1, 2))
#' 
#' stereoplot()
#' angelier(angelier1990$KAM, col = assign_col_d(res$assignment))
#' legend('topleft', title = 'Phases', fill = cols, legend = 1:3)
#' 
#' 
#' stereoplot()
#' angelier(angelier1990$KAM, col = 'grey')
#' for(k in seq_len(res$k_opt)){
#'    res_k <- res$phase_results[[k]]
#'    points(res_k$principal_axes, col = cols[k], pch = 16:18, cex = 2) 
#'  }
#' legend('topright', title = "Principal axes", pch = 16:18, col = 'black', legend = c('S1', 'S2', 'S3'))
#'  
#'  dev.off()
slip_inversion_wissi_polyphase <- function(x, 
                                           weights = NULL,
                                           k_max = 4L,
                                           sigma_K_deg = 30,
                                           seed = NULL, ...){
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  
  res <- wissi_polyphase(normals, slips, 
                  weights = weights,
                  k_max = k_max,
                  sigma_K_deg = sigma_K_deg,
                  seed = seed, ...)
  
  res$phase_results <- lapply(res$phase_results, function(i){
    s <- sigma2stress(i$sigma)
    i$principal_axes <- s$axes
    i$principal_vals <- s$vals
    i$sigma1 <- i$sigma2 <- i$sigma3 <- i$eigenvalues <- i$Phi <- NULL
    
    i$stress_shape <- stress_shape(i$sigma)
    
    i$SHmax <- SH(i$principal_axes[1, ], i$principal_axes[2, ], i$principal_axes[3, ], R = i$stress_shape$R)
    
    i$slips_corrected <- as.Vec3(i$slips_corrected)
    
    i$misfit <- slip_inversion_misfit(i$sigma, x)
    
    return(i)
  })
  
  return(res)
}


# ===
# BOOTSTRAP UNCERTAINTY  (Yamaji-Sato Section 6)
# ===
wissi_bootstrap <- function(normals, slips,
                            weights = NULL,
                            B = 500L,
                            seed = NULL,
                            ...) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(normals)
  res0 <- wissi(normals, slips, weights, ...)
  y0 <- res0$y
  
  resample_one <- function(b) {
    idx <- sample(N, N, replace = TRUE)
    w_b <- if (!is.null(weights)) weights[idx] else NULL
    rb <- tryCatch(
      wissi(
        normals[idx, , drop = FALSE],
        slips[idx, , drop = FALSE],
        w_b, ...
      ),
      error = function(e) NULL
    )
    if (is.null(rb)) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    c(
      angular_stress_distance(y0, rb$y),
      orife_lisle_distance(y0, rb$y),
      michael_distance(y0, rb$y)
    )
  }
  
  boot_mat <- vapply(seq_len(B), resample_one, numeric(3L))
  
  list(
    optimal        = res0,
    thetas         = boot_mat[1, ],
    dispersion = mean(boot_mat[1, ], na.rm = TRUE),
    sd         = sd(boot_mat[1, ], na.rm = TRUE),
    D_bar          = mean(boot_mat[2, ], na.rm = TRUE),
    DM_bar         = mean(boot_mat[3, ], na.rm = TRUE)
  )
}


#' Bootstrap uncertainty for a WISSI result.
#'
#' Yields `n_iter` stress tensors from resampled datasets. The dispersion
#' Theta-bar on S^5 approximates the data noise level (Eq. 37: Theta ~ d-bar).
#'
#' @inheritParams slip_inversion_yamaji_sato_boot
#' @param n_iter     Number of bootstrap replicates. Default `500`.
#' @param seed  Optional RNG seed.
#' @param ...   Additional arguments passed to [slip_inversion_wissi()].
#' 
#' @family wissi 
#'
#' @return A named list with:
#' \describe{
#'   \item{`optimal`}{slip_inversion_wissi() result for the full dataset}
#'   \item{`thetas`}{length-B vector of angular stress distances from optimal}
#'   \item{`dispersion`}{mean Theta (approximates noise level p of data)}
#'   \item{`sd`}{standard deviation of Theta values}
#'   \item{`D_bar`}{mean Orife-Lisle distance from optimal}
#'   \item{`DM_bar`}{mean Michael distance from optimal}
#'   }
#'   
#' @examples
#' res <- slip_inversion_wissi_boot(angelier1990$AVB, n_iter = 4)
#' 
#' stereoplot()
#' angelier(angelier1990$KAM, col = 'grey')
#' lines(res$optimal$principal_axes, res$sd, col = 2:4)
#' points(res$optimal$principal_axes, pch = 16:18, cex = 2, col= 2:4)
#' text(res$optimal$principal_axes, 
#'  label = rownames(res$optimal$principal_axes), col= 2:4, adj = -.25)
slip_inversion_wissi_boot <- function(x, n_iter = 500L, seed = NULL, ...){
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  res <- wissi_bootstrap(normals, slips, B = n_iter, seed = seed, ...)
  
  s <- sigma2stress(res$optimal$sigma)
  res$optimal$principal_axes <- s$axes
  res$optimal$principal_vals <- s$vals
  res$optimal$sigma1 <- res$optimal$sigma2 <- res$optimal$sigma3 <- res$optimal$eigenvalues <- res$optimal$Phi <- NULL
  
  res$optimal$stress_shape <- stress_shape(res$optimal$sigma)
  
  res$optimal$SHmax <- SH(res$optimal$principal_axes[1, ], res$optimal$principal_axes[2, ], res$optimal$principal_axes[3, ], R = res$optimal$stress_shape$R)
  
  res$optimal$slips_corrected <- as.Vec3(res$optimal$slips_corrected)
  
  res$optimal$misfit <- slip_inversion_misfit(res$optimal$sigma, x)
  
  return(res)
}