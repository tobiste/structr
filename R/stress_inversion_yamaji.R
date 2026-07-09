# Fault slip inversion: Yamaji and Sato (2006) method
#
# Reference: Yamaji, A. and Sato, K. (2006). Distances for the solutions of
#   stress tensor inversion in relation to misfit angles that accompany the
#   solutions. Geophys. J. Int., 167, 933-942. doi:10.1111/j.1365-246X.2006.03188.x
#
# Theory summary
# ---
# A reduced stress tensor is represented as a unit 5-vector x on the 5-sphere
# S^5, embedded in the 6D space of symmetric tensors via the y-vector (Eq. 6):
#
#   y = (1/sqrt(2)) * (s11, s22, s33, sqrt(2)*s23, sqrt(2)*s31, sqrt(2)*s12)
#
# with normalisation constraints (Eqs. 11-12):
#   sigmaI  = s11 + s22 + s33 = 0          (deviator)
#   sigmaII = (s11^2+s22^2+s33^2)/2
#             + s23^2+s31^2+s12^2 = 1      (unit norm)
#
# These constraints confine y to the intersection of a hyperplane and the
# 6-sphere, which defines the 5-sphere S^5.
#
# For each fault datum {n, slip}, define the 6-vector eps' (Eq. 8, derived
# from Eq. 4) such that the Wallace-Bott condition becomes:
#
#   eps'_i . y = 0    (shear stress perpendicular to b_i = n_i x slip_i)
#
# The optimal y minimises sum_i (eps'_i . y)^2, i.e. the smallest-eigenvalue
# eigenvector of M = sum_i eps'_i (x) eps'_i.
#
# The computation is carried out in the 5D deviator subspace via an explicit
# orthonormal basis B (6x5), constructed as the complement of the deviator
# null direction h = (1,1,1,0,0,0)/sqrt(3), to avoid degeneracy in the 6x6
# projected matrix.
#
# The sense condition (Eq. 5 / Eq. 10) selects between the two antipodal
# solutions via a majority vote on tau_i . slip_i > 0 (positive work).
#
# Slip convention
# ---
# normals : N x 3 unit fault plane normals
# slips   : N x 3 unit slip direction vectors, parallel to the shear traction
#           (i.e. pointing in the direction of hanging-wall motion). This is
#           OPPOSITE to the paper's v convention (which points against motion).
#
# Distance metrics (Section 3-4)
# ---
# Angular stress distance Theta (Eq. 19): great-circle distance on S^5,
#   Theta = arccos(|y1 . y2|), range [0, 90] deg for reduced tensors.
#
# Orife-Lisle distance D (Eq. 26): Euclidean chord on S^5,
#   D = 2*sin(Theta/2), range [0, 2].
#
# Michael distance DM (Eq. 34): DM = 1 - cos(Theta), range [0, 2].
#
# Relationship to misfit angles (Eq. 37): d-bar approx Theta, where d-bar is
# the mean angular misfit between predicted and observed slip directions.


# Internal helpers

# 6D y-vector from a 3x3 stress tensor (Eq. 6)
.sigma_to_y6 <- function(sigma) {
  c(
    sigma[1, 1], sigma[2, 2], sigma[3, 3],
    sqrt(2) * sigma[2, 3], sqrt(2) * sigma[3, 1], sqrt(2) * sigma[1, 2]
  ) / sqrt(2)
}

# 3x3 stress tensor from a 6D y-vector (inverse of Eq. 6)
.y6_to_sigma <- function(y) {
  sv <- y * sqrt(2)
  sigma <- matrix(0, 3, 3)
  diag(sigma) <- sv[1:3]
  sigma[2, 3] <- sigma[3, 2] <- sv[4] / sqrt(2)
  sigma[3, 1] <- sigma[1, 3] <- sv[5] / sqrt(2)
  sigma[1, 2] <- sigma[2, 1] <- sv[6] / sqrt(2)
  sigma
}

# eps'_i : 6D datum vector for the Wallace-Bott constraint (derived from Eq. 4)
# b_i = n_i x slip_i (unit vector in fault plane, perpendicular to slip)
# Satisfies: eps'_i . y = sum_{jk} sigma_jk * (b_j*n_k + b_k*n_j) / 2 = 0
.eps6_prime <- function(slip, n) {
  b <- c(
    n[2] * slip[3] - n[3] * slip[2],
    n[3] * slip[1] - n[1] * slip[3],
    n[1] * slip[2] - n[2] * slip[1]
  )
  c(
    b[1] * n[1] * sqrt(2), b[2] * n[2] * sqrt(2), b[3] * n[3] * sqrt(2),
    b[2] * n[3] + b[3] * n[2],
    b[1] * n[3] + b[3] * n[1],
    b[1] * n[2] + b[2] * n[1]
  )
}

# Orthonormal basis B for the 5D deviator subspace (complement of h)
# h = (1,1,1,0,0,0)/sqrt(3) is the deviator null direction (Eq. 15)
# B is a 6x5 matrix; computed once at load time
.make_basis_B <- function() {
  h <- c(1, 1, 1, 0, 0, 0) / sqrt(3)
  Q_full <- qr.Q(qr(cbind(h, diag(6)))) # 6x6 orthonormal, first col = h
  Q_full[, 2:6] # 6x5 complement basis
}
.B <- .make_basis_B()

# Sense resolution: choose y or -y such that tau_i . slip_i > 0 for majority
# (positive work condition, Eq. 10 in paper's convention)
.resolve_sense <- function(y_raw, normals, slips) {
  sigma <- .y6_to_sigma(y_raw)
  dots <- vapply(seq_len(nrow(normals)), function(i) {
    tau <- sigma %*% normals[i, ]
    tau <- tau - sum(tau * normals[i, ]) * normals[i, ]
    sum(tau * slips[i, ])
  }, numeric(1L))
  if (sum(dots < 0) > sum(dots > 0)) -y_raw else y_raw
}


yamaji_sato <- function(normals, slips, wt) {
  N <- nrow(normals)

  # Input validation
  if (N < 4L) {
    stop("At least 4 fault slip measurements are required.")
  }
  stopifnot(N == length(wt))

  # --- Linear step: eigenvector problem in 5D (Section 4.2.7 of Pascal 2022) ---
  #
  # Accumulate M5 = sum_i w_i * (B^T eps'_i)(B^T eps'_i)^T   (5x5)
  # The smallest-eigenvalue eigenvector x5 of M5 satisfies
  # eps'_i . (B x5) ~ 0 for all i, i.e. the optimal 5D stress vector.
  M5 <- matrix(0, 5, 5)
  for (i in seq_len(N)) {
    ep5 <- drop(crossprod(.B, .eps6_prime(slips[i, ], normals[i, ]))) # 5-vector
    M5 <- M5 + wt[i] * tcrossprod(ep5)
  }

  eig <- eigen(M5, symmetric = TRUE) # eigenvalues in decreasing order
  x5 <- eig$vectors[, 5L] # smallest-eigenvalue eigenvector

  # Map back to 6-space and enforce unit norm
  y_raw <- drop(.B %*% x5)
  y_raw <- y_raw / sqrt(sum(y_raw^2))

  # --- Nonlinear step: sense resolution (Eq. 10 / Eq. 5) ---
  .resolve_sense(y_raw, normals, slips)
}

#' Stress tensor inversion via the Yamaji and Sato (2006) eigenvector method.
#'
#' @param x object of class `"Pair"` or `"Fault"` with at least 4 rows.
#' @inheritParams slip_inversion_angelier
#'
#' @note Note on opposite-tensor ambiguity
#'
#' The Wallace-Bott condition eps'_i . y = 0 is satisfied by both y and -y,
#'  corresponding to stress tensors.
#'  The sense condition (positive work) selects the physically meaningful sign,
#'  but requires that the majority of slip vectors have correct sense. If more
#'  than half the slips are recorded in the wrong sense, the opposite tensor
#'  will be returned. Check alpha_signed_deg and suspected_flipped accordingly.
#'
#'  @references
#'  Yamaji, A., & Sato, K. (2006). Distances for the solutions of stress tensor
#' inversion in relation to misfit angles that accompany the solutions.
#' Geophysical Journal International, 167(2), 933–942.
#' https://doi.org/10.1111/j.1365-246X.2006.03188.x
#'
#' @returns Same output as [slip_inversion()] plus
#' \describe{
#' \item{y}{6D unit y-vector on S^5 representing the tensor}
#' \item{alpha}{per-fault angular misfit (unsigned, 0-90&deg;)}
#' \item{mean_alpha}{mean angular misfit across all faults}
#' }
#' @export
#'
#' @family stress-inversion
#'
#' @examples
#' nx <- length(angelier1990)
#' par(mfrow = c(1, nx))
#'
#' invisible(lapply(seq_len(nx), function(i) {
#'   # inversion
#'   x <- angelier1990[[i]]
#'   res <- slip_inversion_yamaji_sato(x)
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
slip_inversion_yamaji_sato <- function(x, weights = NULL, flip = FALSE, friction = 0.6) {
  tsign <- if (flip) -1 else 1

  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()
  N <- nrow(normals)

  # Weights: normalise and impute NAs with observed mean
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

  y_opt <- yamaji_sato(normals, slips, wt)

  # --- Post-processing ---
  TR <- .y6_to_sigma(y_opt) * tsign
  # eig_s  <- eigen(sigma, symmetric = TRUE)
  # eigval <- eig_s$values                  # sigma1 >= sigma2 >= sigma3
  # eigvec <- eig_s$vectors
  # Phi    <- (eigval[1] - eigval[2]) / (eigval[1] - eigval[3])

  p <- tau2stress(TR)
  shape <- stress_shape(TR)
  misfit <- slip_inversion_misfit(TR, x)

  # Per-fault angular misfit (unsigned: 0-90 deg, the standard line metric)
  alpha_deg <- vapply(seq_len(N), function(i) {
    tau <- TR %*% normals[i, ]
    tau <- tau - sum(tau * normals[i, ]) * normals[i, ]
    tn <- sqrt(sum(tau^2))
    if (tn < 1e-12) {
      return(NA_real_)
    }
    acos(pmax(-1, pmin(1, abs(sum(tau / tn * slips[i, ]))))) * 180 / pi
  }, numeric(1L))


  # Angle between slip planes and sigma 1
  theta <- vapply(seq_len(N), function(i) {
    angle(Plane(x[i, ]), p$principal_axes[1, ])
  }, numeric(1))


  # Per-fault signed misfit (0-180 deg): values > 90 indicate probable
  # slip sense recording error
  # alpha_signed <- vapply(seq_len(N), function(i) {
  #   tau <- sigma %*% normals[i, ]
  #   tau <- tau - sum(tau * normals[i, ]) * normals[i, ]
  #   tn  <- sqrt(sum(tau^2))
  #   if (tn < 1e-12) return(NA_real_)
  #   acos(pmax(-1, pmin(1, sum(tau / tn * slips[i, ])))) * 180 / pi
  # }, numeric(1L))


  # Theoretically resolved shear stress on plane
  sigma_s_mean <- mean(abs(shear_stress(p$sigma_vals[1], p$sigma_vals[3], theta)))

  shearnorm <- tau2shearnorm(TR, x, friction = friction)
  tendency <- tau2tendency(TR, x, friction = friction)

  # sigma_s_mean <- mean(abs(shearnorm))

  SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(TR, symmetric = TRUE)$vectors),
    error = function(e) {
      SH(p$principal_axes[1, ], p$principal_axes[2, ], p$principal_axes[3, ], R = shape$R)
    }
  )

  pfaults <- principal_fault(p$principal_axes[1, ], p$principal_axes[3, ], friction)

  list(
    stress_tensor = TR,
    # tensor_params = tensor_params,
    principal_axes = p$principal_axes,
    principal_vals = p$sigma_vals,
    principal_faults = pfaults,
    stress_shape = shape,
    tau_mean = sigma_s_mean,
    stress_components = cbind(shearnorm, tendency),
    misfit = misfit,
    SHmax = SHmax,
    alpha = alpha_deg,
    # alpha_signed_deg = alpha_signed,
    mean_alpha = mean(alpha_deg, na.rm = TRUE),
    y = y_opt,
    method = "yamaji_sato"
  )
}


# Distance metrics (Section 3-4 of Yamaji and Sato 2006)
# All functions accept 6D y-vectors as returned by yamaji_sato().

# Angular stress distance Theta (Eq. 19): great-circle distance on S^5.
# Range [0, 90] deg. The abs() handles the opposite-tensor ambiguity
# (y and -y represent equivalent stress states for the purposes of distance).
# Eq. 37: Theta approximately equals the mean misfit angle d-bar.
angular_stress_distance <- function(y1, y2) {
  acos(pmax(-1, pmin(1, abs(sum(y1 * y2))))) * 180 / pi
}

# Orife-Lisle distance D (Eq. 26): Euclidean chord length on S^5.
# D = 2*sin(Theta/2), range [0, sqrt(2)] for reduced tensors.
# Equivalent to |x1 - x2| on the 5-sphere.
orife_lisle_distance <- function(y1, y2) {
  theta <- angular_stress_distance(y1, y2) * pi / 180
  2 * sin(theta / 2)
}

# Michael distance DM (Eq. 34): DM = 1 - cos(Theta), range [0, 1].
# DM = D^2 / 2 (Eq. 35). Less sensitive than Theta or D for small perturbations.
michael_distance <- function(y1, y2) {
  theta <- angular_stress_distance(y1, y2) * pi / 180
  1 - cos(theta)
}


#' Uncertainties of direction stress inversion after Yamaji and Sato (2006)
#'
#' Bootstrap resampling to evaluate solution precision (Section 6).
#' Yields B stress tensors from resampled datasets. The dispersion of these
#' tensors on S^5 approximates the noise level of the data (Eq. 37).
#'
#' @inheritParams slip_inversion_yamaji_sato
#' @param n_boot integer. Number of bootstrap replicates
#'
#' @returns list:
#' \describe{
#' \item{optimal}{[slip_inversion_yamaji_sato()] result for the full dataset}
#' \item{thetas}{length-B vector of angular stress distances from optimal}
#' \item{dispersion_deg}{ mean angular stress distance (Theta-bar); approximates
#' the noise level p of the data (Fig. 8 of paper)}
#' \item{sd_deg}{standard deviation of Theta values}
#' \item{D_bar}{mean Orife-Lisle distance from optimal}
#' \item{DM_bar}{mean Michael distance from optimal}
#' }
#' @export
#'
#' @seealso [slip_inversion_yamaji_sato()]
#'
#' @examples
#' slip_inversion_yamaji_sato_boot(angelier1990$AVB)
slip_inversion_yamaji_sato_boot <- function(x, weights = NULL,
                                            n_boot = 500L) {
  # if (!is.null(seed)) set.seed(seed)

  normals <- Vec3(Plane(x)) |> unclass()
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips) |> unclass()

  N <- nrow(normals)

  # Weights: normalise and impute NAs with observed mean
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

  res0 <- slip_inversion_yamaji_sato(x, wt)
  y0 <- res0$y

  thetas <- vapply(seq_len(n_boot), function(b) {
    idx <- sample(N, N, replace = TRUE)
    y <- tryCatch(
      yamaji_sato(normals[idx, ], slips[idx, ], wt[idx]),
      error = function(e) NULL
    )
    if (is.null(y)) {
      return(NA_real_)
    }
    angular_stress_distance(y0, y)
  }, numeric(1L))

  D_vals <- orife_lisle_distance(y0, y0) # placeholder; computed per resample below
  DM_vals <- michael_distance(y0, y0)

  # Recompute per-replicate distances properly
  all_y <- lapply(seq_len(n_boot), function(b) {
    idx <- sample(N, N, replace = TRUE)
    y <- tryCatch(
      yamaji_sato(normals[idx, ], slips[idx, ], wt[idx]),
      error = function(e) NULL
    )
    if (is.null(y)) {
      return(NULL)
    } else {
      y
    }
  })

  valid <- !sapply(all_y, is.null)
  thetas2 <- sapply(all_y[valid], function(y) angular_stress_distance(y0, y))
  D_bar <- mean(sapply(all_y[valid], function(y) orife_lisle_distance(y0, y)))
  DM_bar <- mean(sapply(all_y[valid], function(y) michael_distance(y0, y)))


  list(
    optimal        = res0,
    thetas         = thetas2,
    dispersion_deg = mean(thetas2, na.rm = TRUE),
    sd_deg         = sd(thetas2, na.rm = TRUE),
    D_bar          = D_bar,
    DM_bar         = DM_bar
  )
}


# Utility: build normalised reduced stress tensor from Phi and principal axes.
# R_mat: 3x3, columns = (sigma1_dir, sigma2_dir, sigma3_dir) unit vectors.
# Follows Eq. 13-14 of Yamaji and Sato (2006).
build_reduced_sigma <- function(Phi, R_mat) {
  Rsc <- sqrt(3 * Phi^2 - 3 * Phi + 3)
  R_mat %*% diag(c(2 - Phi, 2 * Phi - 1, -Phi - 1) / Rsc) %*% t(R_mat)
}
