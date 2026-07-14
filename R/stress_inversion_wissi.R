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
signal_to_weights <- function(signal,
                              scale_fn = function(x) 1 / x^2,
                              pessimistic = FALSE) {
  stopifnot(is.numeric(signal), any(!is.na(signal)))
  fill <- if (pessimistic) {
    quantile(signal, 0.75, na.rm = TRUE)
  } else {
    mean(signal, na.rm = TRUE)
  }
  complete <- ifelse(is.na(signal), fill, signal)
  scale_fn(complete)
}

#' Combine multiple weight vectors multiplicatively and normalise.
#'
#' Each element of ... is a numeric vector of length N.
#' Returns a normalised weight vector (mean = 1).
combine_weights <- function(...) {
  wlist <- list(...)
  if (length(wlist) == 0L) stop("Provide at least one weight vector.")
  lens <- lengths(wlist)
  if (length(unique(lens)) != 1L) {
    stop("All weight vectors must have the same length.")
  }
  w <- Reduce(`*`, wlist)
  if (any(w < 0)) stop("Combined weights must be non-negative.")
  if (!any(w > 0)) stop("At least one weight must be positive.")
  w / mean(w)
}


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
#   2. Build the Gaussian affinity matrix K_ij = exp(-D_ij^2 / 2*sigma_K^2).
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
  km <- kmeans(U_norm, centers = k_opt, nstart = 25L, iter.max = 200L)
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
#' Combines the four classic fault slip inversion algorithms into a single
#' coherent framework operating in the Yamaji-Sato 5-sphere sigma-space.
#'
#' @param normals         N x 3 matrix of unit fault plane normals.
#' @param slips           N x 3 matrix of unit slip direction vectors
#'                        (parallel to shear traction; direction of
#'                        hanging-wall motion).
#' @param weights         Optional length-N weight vector. NAs imputed with
#'                        observed mean. Normalised internally (mean = 1).
#'                        Use signal_to_weights() and combine_weights() to
#'                        build from field ranks, measurement errors, and
#'                        prior RUP values.
#' @param sigma_alpha_deg Estimated slip direction measurement error in
#'                        degrees. Used in Stage 4 analytic uncertainty.
#'                        Default 10.
#' @param gamma_max       Maximum sense annealing sharpness parameter.
#'                        Higher values commit more strongly to the predicted
#'                        slip sense. Default 10. Set to 0 to disable sense
#'                        annealing (fully sense-agnostic, like Hansen 2013).
#' @param n_anneal        Number of annealing steps (outer loop of Stage 3).
#'                        Default 8.
#' @param max_iter        Maximum inner iterations per annealing step. Default 50.
#' @param tol_deg         Convergence tolerance in angular stress distance
#'                        (degrees). Default 1e-4.
#' @param run_stage4      Logical. Compute analytic uncertainty (Stage 4).
#'                        Default TRUE.
#'
#' @return A named list with:
#'   sigma             : 3x3 reduced stress tensor (Cartesian frame)
#'   y                 : 6D unit y-vector on S^5
#'   sigma1/2/3        : unit vectors of principal stress axes (max to min)
#'   eigenvalues       : eigenvalues of sigma (decreasing)
#'   Phi               : shape ratio (sigma1-sigma2)/(sigma1-sigma3) in \eqn{[0,1]}
#'   alpha_deg         : per-fault angular misfit (unsigned, 0-90 deg)
#'   alpha_signed_deg  : per-fault signed misfit (0-180 deg)
#'   mean_alpha        : mean angular misfit across all faults (deg)
#'   suspected_flipped : row indices where alpha_signed > 90 deg
#'   n_flipped_sense   : number of faults whose sense was corrected in Stage 3
#'   slips_corrected   : sense-corrected slip matrix used in final inversion
#'   mu                : per-fault magnitude weights from Stage 2/3
#'   phi_sense         : per-fault tanh sense confidence from Stage 3
#'   eigenvalue_gap    : lambda_2 - lambda_1 of M5 (condition number proxy)
#'   M5_eigvals        : all 5 eigenvalues of final M5
#'   unc               : Stage 4 uncertainty list (if run_stage4 = TRUE):
#'     Cov5            : 5x5 covariance matrix in sigma-space
#'     Cov_y6          : 6x6 covariance matrix (y-space)
#'     eigval_gap      : eigenvalue gap (same as above)
#'     cov_eigvals     : eigenvalues of Cov5
#'     sigma1_unc_deg  : approx 1-sigma uncertainty on sigma1 orientation
#'     Phi_unc         : approx 1-sigma uncertainty on Phi
#'   n_iter_total      : total number of inner iterations
wissi <- function(normals,
                  slips,
                  weights = NULL,
                  sigma_alpha_deg = 10,
                  gamma_max = 10,
                  n_anneal = 8L,
                  max_iter = 50L,
                  tol_deg = 1e-4,
                  run_stage4 = TRUE) {
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

  # --- Stage 1: Weighted eigenvector initialisation ---
  s1 <- .wissi_stage1(normals, slips, wt)

  # --- Stage 2: Magnitude reweighting (Mostafa in sigma-space) ---
  s2 <- .wissi_stage2(normals, slips, wt, s1$y,
    max_iter = max_iter, tol_deg = tol_deg
  )

  # --- Stage 3: Sense annealing ---
  s3 <- .wissi_stage3(normals, slips, wt, s2$y,
    gamma_max = gamma_max,
    n_anneal  = n_anneal,
    max_iter  = max_iter,
    tol_deg   = tol_deg
  )

  y_opt <- s3$y
  slips_corr <- s3$slips_corrected
  wt_final <- s3$wt_final

  # --- Stage 4: Analytic uncertainty ---
  unc <- if (run_stage4) {
    .wissi_stage4(normals, slips_corr, wt_final, y_opt,
      sigma_alpha_deg = sigma_alpha_deg
    )
  } else {
    NULL
  }

  # --- Post-processing ---
  dec <- .decompose_y(y_opt)
  sigma <- dec$sigma

  principal_axes <- rbind(dec$sigma1, dec$sigma2, dec$sigma3)
  principal_vals <- dec$eigenvalues
  rownames(principal_axes) <- names(principal_vals) <- c("sigma1", "sigma2", "sigma3")
  
  a_deg <- .alpha_deg(sigma, normals, slips_corr)
  a_signed <- .alpha_signed_deg(sigma, normals, slips_corr)

  list(
    stress_tensor             = sigma,
    y                 = y_opt,
    principal_axes    = principal_axes,
    principal_vals    = principal_vals,
    #Phi               = dec$Phi,
    alpha             = a_deg,
    alpha_signed      = a_signed,
    mean_alpha        = mean(a_deg, na.rm = TRUE),
    suspected_flipped = which(a_signed > 90),
    n_flipped_sense   = s3$n_flipped,
    slips_corrected   = slips_corr,
    mu                = s3$mu,
    phi_sense         = s3$phi,
    eigenvalue_gap    = s3$eigenvalue_gap,
    M5_eigvals        = s3$eigvals,
    unc               = unc,
    n_iter_total      = s3$n_iter_total
  )
}

#' @title Weighted Iterative Sigma-Space Inversion (WISSI)
#'
#' @description Combines the four classic fault slip inversion algorithms into a single
#' coherent framework operating in the Yamaji-Sato 5-sphere sigma-space.
#'
#' @inheritParams slip_inversion_yamaji_sato
#' @param sigma_alpha_deg Estimated slip direction measurement error in
#'                        degrees. Used in Stage 4 analytic uncertainty.
#'                        Default `10.`
#' @param gamma_max       Maximum sense annealing sharpness parameter.
#'                        Higher values commit more strongly to the predicted
#'                        slip sense. Default `10.` Set to `0` to disable sense
#'                        annealing (fully sense-agnostic, like Hansen 2013).
#' @param n_anneal        Number of annealing steps (outer loop of Stage 3).
#'                        Default `8`.
#' @param max_iter        Maximum inner iterations per annealing step. Default `50.`
#' @param tol_deg         Convergence tolerance in angular stress distance
#'                        (degrees). Default `1e-4`.
#' @param run_stage4      Logical. Compute analytic uncertainty (Stage 4).
#'                        Default `TRUE.`
#'
#' @returns A named list with: \describe{
#'   \item{`stress_tensor`}{3x3 reduced stress tensor (Cartesian frame)}
#'   \item{`y`}{6D unit y-vector on \eqn{S^5}}
#'   \item{`principal_axes`}{unit vectors of principal stress axes (max to min)}
#'   \item{`principal_vals`}{eigenvalues of `stress_tensor` (decreasing)}
#'   \item{`alpha`}{per-fault angular misfit (unsigned, 0-90&deg;)}
#'   \item{`alpha_signed`}{per-fault signed misfit (0-180&deg;)}
#'   \item{`mean_alpha`}{mean angular misfit across all faults (&deg;)}
#'   \item{`suspected_flipped`}{row indices where `alpha_signed` > 90&deg;}
#'   \item{`n_flipped_sense`}{number of faults whose sense was corrected in Stage 3}
#'   \item{`slips_corrected`}{sense-corrected slip matrix used in final inversion}
#'   \item{`mu`}{per-fault magnitude weights from Stage 2/3}
#'   \item{`phi_sense`}{per-fault tanh sense confidence from Stage 3}
#'   \item{`eigenvalue_gap`}{\eqn{\lambda_2 - \lambda_1} of \eqn{M5} (condition number proxy)}
#'   \item{`M5_eigvals`}{all 5 eigenvalues of final \eqn{M5}}
#'   \item{`unc`}{Stage 4 uncertainty list (if `run_stage4 = TRUE`): 
#'     `Cov5`: 5x5 covariance matrix in sigma-space;
#'     `Cov_y6`: 6x6 covariance matrix (`y`-space);
#'     `eigval_gap`: eigenvalue gap (same as above);
#'     `cov_eigvals`: eigenvalues of `Cov5`;
#'     `sigma1_unc`: approx 1\eqn{\sigma} uncertainty on \eqn{\sigma_1} orientation
#'     `Phi_unc`: approx 1\eqn{sigma} uncertainty on \eqn{\phi}}
#'   \item{`n_iter_total`}{total number of inner iterations}
#' }
#'   
#' @family stress-inversion
#' 
#' @references Stephan (in prep.)
#'   
#' @export
#'
#' @examples
#' nx <- length(angelier1990)
#' par(mfrow = c(1, nx))
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
slip_inversion_wissi <- function(x, weights = NULL,
                                 sigma_alpha_deg = 10,
                                 gamma_max = 10,
                                 n_anneal = 8L,
                                 max_iter = 50L,
                                 tol_deg = 1e-4,
                                 run_stage4 = TRUE) {
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  res <- wissi(normals, slips, weights = weights, sigma_alpha_deg = sigma_alpha_deg,
               gamma_max = gamma_max, n_anneal = n_anneal, max_iter = max_iter,
               tol_deg = tol_deg, run_stage4 = run_stage4)
  
  res$principal_axes <- as.Vec3(res$principal_axes) |> Line()
  res$stress_shape <- stress_shape(res$stress_tensor)
  
  res$SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(res$stress_tensor, symmetric = TRUE)$vectors),
    error = function(e) {
      SH(res$principal_axes[1, ], res$principal_axes[2, ], res$principal_axes[3, ], R = res$stress_shape$R)
    }
  )
  
  res$slips_corrected <- as.Vec3(res$slips_corrected)
  
  res$misfit <- slip_inversion_misfit(res$sigma, x)
  
  return(res)
}


# WISSI_POLYPHASE — Polyphase separation wrapper

#
#' Polyphase stress inversion via spectral clustering on \eqn{S^5} (Stage 5).
#'
#' Identifies k stress phases automatically using the eigengap heuristic,
#' then runs wissi() on each phase subset.
#'
#' @inheritParams slip_inversion_wissi
#' @param k_max       Maximum number of phases to consider. Default 4.
#' @param sigma_K_deg Affinity bandwidth in degrees of angular stress
#'                    distance. Faults within this distance are considered
#'                    similar. Default 30. Increase to merge nearby phases,
#'                    decrease to split them.
#' @param seed        Optional RNG seed for k-means reproducibility.
#' @param ...         Additional arguments passed to wissi() for each phase.
#'
#' @return A named list with:
#'   assignment    : integer vector of length N (phase label 1..k per fault)
#'   k_opt         : number of phases identified
#'   gaps          : Laplacian eigenvalue gaps (eigengap criterion)
#'   phase_results : list of k wissi() results, one per phase
#'   D_mat         : N x N pairwise ASD matrix between fault poles (degrees)
#' @export
#' @examples
#' res <- wissi_polyphase(angelier1990$TYM)
#' res$k_opt
wissi_polyphase <- function(normals, slips,
                            weights = NULL,
                            k_max = 4L,
                            sigma_K_deg = 30,
                            seed = NULL,
                            ...) {
  
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  
  # if (!is.matrix(normals) || !is.matrix(slips)) {
  #   stop("'normals' and 'slips' must be matrices.")
  # }
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


# ===
# BOOTSTRAP UNCERTAINTY  (Yamaji-Sato Section 6)
# ===
#
#' Bootstrap uncertainty for a wissi() result.
#'
#' Yields B stress tensors from resampled datasets. The dispersion
#' Theta-bar on S^5 approximates the data noise level (Eq. 37: Theta ~ d-bar).
#'
#' @param normals,slips,weights As for wissi().
#' @param B     Number of bootstrap replicates. Default 500.
#' @param seed  Optional RNG seed.
#' @param ...   Additional arguments passed to wissi().
#'
#' @return A named list with:
#'   optimal        : wissi() result for the full dataset
#'   thetas         : length-B vector of angular stress distances from optimal
#'   dispersion : mean Theta (approximates noise level p of data)
#'   sd         : standard deviation of Theta values
#'   D_bar          : mean Orife-Lisle distance from optimal
#'   DM_bar         : mean Michael distance from optimal
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
