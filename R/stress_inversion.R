#' Stress Inversion for Fault-Slip Data
#' 
#' Convenient function for linear and direction inversion of fault-slip data to 
#' derive the reduced stress tensor.
#'
#' @param x `"Fault"` object where the rows are the observations, and the columns the coordinates.
#' @param method character. The inversion algorithm. One of `"michael"` (the default) for a bootstrapped linear inversion, 
#' or `"angelier"` for a iterative direct inversion.
#' @param ... arguments passed to [slip_inversion_angelier()] or [slip_inversion_michael()]
#'
#' @returns list. See [slip_inversion_angelier()] and [slip_inversion_michael()]
#'
#' @family stress-inversion
#' 
#' @references 
#' Angelier, J. (1990). Inversion of field data in fault tectonics to obtain 
#' the regional stress—III. A new rapid direct inversion method by analytical 
#' means. Geophys. J. Int, 103, 363–376. <https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>
#' 
#' Michael, A. J. (1984). Determination of stress from slip data: Faults and 
#' folds. Journal of Geophysical Research: Solid Earth, 89(B13), 11517–11526. 
#' <https://doi.org/10.1029/JB089iB13p11517>
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#' # Use Angelier examples
#' par(mfrow = c(1, length(angelier1990)))
#' invisible(lapply(angelier1990, function(x){
#' 
#' # Linear inversion
#' TYM_michael <- slip_inversion(x, method = "michael")
#' 
#' # Direct inversion
#' TYM_angelier <- slip_inversion(x, method = "angelier")
#' 
#' stereoplot(guides = FALSE)
#' fault_plot(x, col = "gray80")
#' points(TYM_michael$principal_axes, pch = 3, col = 1:3)
#' points(TYM_angelier$principal_axes, pch = 2, col = 1:3)
#' legend("topleft", col = 1, legend = c("Michael (1984)", "Angelier (1990)"), pch = c(3, 2))
#' }))
slip_inversion <- function(x, method = c("michael", "angelier"), ...) {
  method <- match.arg(method)

  if (method == "angelier") {
    slip_inversion_angelier(x, ...)
  } else {
    slip_inversion_michael(x, ...)
  }
}

# Michael (1894) ---------------------------------------------------------------

#' Stress Inversion for Fault-Slip Data after Michael (1984)
#'
#' Linear stress inversion (based on Michael, 1984) determines the orientation
#' of the principal stresses from fault slip data.
#' Confidence intervals are estimated by bootstrapping.
#' This inversion is simplified by the assumption that the magnitude of the
#' tangential traction on the various fault planes, at the time of rupture, is similar.
#'
#' @inheritParams slip_inversion
#' @param n_iter integer. Number of bootstrap samples (10 by default)
#' @param conf.level numeric. Confidence level of the interval (0.95 by default)
#' @param friction numeric. Coefficient of friction (0.6 by default)
#' @param ... optional parameters passed to [confidence_ellipse()]
#'
#' @returns list
#'  \describe{
#'  \item{`stress_tensor`}{matrix. Best-fit devitoric stress tensor}
#'  \item{`principal_axes`}{`"Line"` objects. Orientation of the principal stress axes}
#'  \item{`principal_axes_conf`}{list containing the confidesnce ellipses for the 3 principal stress vectors. See [confidence_ellipse()] for details.}
#'  \item{`principal_vals`}{numeric. The proportional magnitudes of the principal stress axes given by the eigenvalues of the stress tensor: \eqn{\sigma_1}, \eqn{\sigma_2}, and \eqn{\sigma_3}}
#'  \item{`principal_vals_conf`}{3-column vector containing the lower and upper margins of the confidence interval of the principal vals}
#'  \item{`principal_fault`}{Principal fault planes as `"Fault"` objects.}
#'  \item{`SHmax`}{numeric. Direction of maximum horizontal stress (in degrees)}
#'  \item{`SHmax_CI`}{numeric. Confidence interval of `SHmax` angle}
#'  \item{`R`}{numeric. Stress shape ratio after Gephart & Forsyth (1984): \eqn{R = (\sigma_1 - \sigma_2)/(\sigma_1 - \sigma_3)}. Values ranging from 0 to 1, with 0 being
#' \eqn{\sigma_1 = \sigma_2} and 1 being \eqn{\sigma_2 = \sigma_3}.}
#'  \item{`R_conf`}{Confidence interval for `R`}
#'  \item{`phi`}{numeric. Stress shape ratio after Angelier (1979): \eqn{\Phi = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)}. Values range between 0 (\eqn{\sigma_2 = \sigma_3}) and 1 (\eqn{\sigma_2 = \sigma_1}).}
#'  \item{`phi_conf`}{Confidence interval for `phi`}
#'  \item{`bott`}{numeric. Stress shape ratio after Bott (1959): \eqn{\R = (\sigma_3 - \sigma_1)/(\sigma_2 - \sigma_1)}. Values range between \eqn{-\infty} and \eqn{+\infty}.}
#'  \item{`bott_conf`}{Confidence interval for `bott`}
#'  \item{`beta`}{numeric. Average angle between the tangential traction predicted by the best stress tensor  and the slip vector on each plane. Should be close to 0.}
#'  \item{`beta_CI`}{numeric. Confidence interval of `beta` angle}
#'  \item{`sigma_s`}{numeric. Average resolved shear stress on each plane. Should be close to 1.}
#'  \item{`fault_data`}{`data.frame` containing the beta angles, the angles between sigma 1 and the plane normal,
#'  the resolved shear and normal stresses, and the slip and dilation tendency on each plane.}
#'  }
#'
#' @references Michael, A. J. (1984). Determination of stress from slip data:
#' Faults and folds. Journal of Geophysical Research: Solid Earth, 89(B13), 11517–11526.
#' \doi{10.1029/JB089iB13p11517}
#'
#' @details
#' The goal of slip inversion is to find the single uniform stress tensor that
#' most likely caused the faulting events. With only slip data to constrain the
#' stress tensor the isotropic component can not be determined, unless
#' assumptions about the fracture criterion are made. Hence inversion will be
#' for the deviatoric stress tensor only.
#' A single fault can not completely constrain the deviatoric stress tensor a,
#' therefore it is necessary to simultaneously solve for a number of faults, so
#' that a single a that best satisfies all of the faults is found.
#'
#' @export
#' @family stress-inversion
#'
#' @importFrom stats t.test
#' @importFrom future.apply future_lapply
#'
#' @seealso [Fault_PT()] for a simple P-T stress analysis,
#'  [SH()] and [SH_from_tensor()] to calculate the azimuth of the maximum horizontal stress;
#'  [Mohr_plot()] for graphical representation of the deviatoric stress tensor.
#'
#' @examples
#' # Use Angelier examples:
#' par(mfrow = c(1, length(angelier1990)))
#' 
#' invisible(lapply(angelier1990, function(x){
#' 
#' # inversion
#' res <- slip_inversion_michael(x, n_iter = 100, n = 1000, res = 100)
#' 
#' # some stress shape
#' R_val <- round(res$R, 2)
#' 
#' # misfit
#' rup_val <- round(res$mean_rup, 2)
#'
#' # Plot the faults (color-coded by RUP%) and show the principal stress axes
#' stereoplot(title = "Bootstrapped linear inversion", guides = FALSE)
#' stereo_shmax(res$SHmax)
#' fault_plot(x, col = assign_col(res$fault_data$rup))
#' stereo_confidence(res$principal_axes_conf$sigma1, col = 2)
#' stereo_confidence(res$principal_axes_conf$sigma2, col = 3)
#' stereo_confidence(res$principal_axes_conf$sigma3, col = 4)
#' text(res$principal_axes, label = rownames(res$principal_axes), col = 2:4, adj = -.25)
#' legend("topleft", col = 2:4, legend = rownames(res$principal_axes), pch = 16)
#' title(sub = bquote(R == .(R_val) ~ "|" ~ bar("RUP") == .(rup_val) * '%'))
#' }))
slip_inversion_michael <- function(x, n_iter = 100L, conf.level = 0.95, friction = 0.6, ...) {
  best.fit <- .slip_inversion_michael(x, friction)
  fault_df <- best.fit$fault_data
  nx <- nrow(x)

  normals <- Fault_plane(x) |>
    Vec3() |>
    unclass()
  slips <- Fault_slip(x) |>
    Vec3() |>
    unclass()

  # bootstrap results
  tau_boot <- future.apply::future_lapply(seq_len(n_iter), function(i) {
    idx <- sample.int(nx, replace = TRUE)
    x_sample <- x[idx, ]
    linear_stress_inversion(normals[idx, ], slips[idx, ])
  }, future.seed = TRUE)

  princ_boot <- lapply(tau_boot, tau2stress)

  # calculate confidence intervals from bootstrap results
  sigma_vec1 <- do.call(rbind, lapply(princ_boot, function(x) {
    x$principal_axes[1, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)

  sigma_vec2 <- do.call(rbind, lapply(princ_boot, function(x) {
    x$principal_axes[2, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)

  sigma_vec3 <- do.call(rbind, lapply(princ_boot, function(x) {
    x$principal_axes[3, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)


  params_boot <- lapply(tau_boot, tau_params)

  R_boot <- vapply(params_boot, function(x) {
    x$R
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  phi_boot <- vapply(params_boot, function(x) {
    x$phi
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  bott_boot <- vapply(params_boot, function(x) {
    x$bott
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  sigma_boot0 <- vapply(princ_boot, function(x) {
    x$sigma_vals
  }, FUN.VALUE = numeric(3)) |> t()
  sigma_boot <- sapply(1:3, function(col) {
    sigma_boot_col <- stats::t.test(sigma_boot0[, col], conf.level = conf.level)
    sigma_boot_col$conf.int
  })
  colnames(sigma_boot) <- names(best.fit$principal_vals)

  beta_CI <- tectonicr::confidence_interval(fault_df$beta, conf.level = conf.level, axial = FALSE)


  SHmax_CI <- vapply(tau_boot, function(x) {
    tryCatch(
      expr = SH_from_tensor(x),
      error = function(e) {
        phi <- tau_params(x)$phi
        principal_axes <- tau2stress(x)$principal_axes
        SH(principal_axes[1, ], principal_axes[2, ], principal_axes[3, ], R = phi)
      }
    )
  }, FUN.VALUE = numeric(1)) |>
    tectonicr::confidence_interval(conf.level = conf.level, axial = TRUE)

  list(
    stress_tensor = best.fit$stress_tensor,
    principal_axes = best.fit$principal_axes,
    principal_axes_conf = list(sigma1 = sigma_vec1, sigma2 = sigma_vec2, sigma3 = sigma_vec3),
    principal_vals = best.fit$principal_vals,
    principal_vals_conf = sigma_boot,
    principal_faults = best.fit$principal_faults,
    SHmax = best.fit$SHmax,
    SHmax_CI = SHmax_CI$conf.interval,
    R = best.fit$R,
    R_conf = R_boot$conf.int,
    phi = best.fit$phi,
    phi_conf = phi_boot$conf.int,
    bott = best.fit$bott,
    bott_conf = bott_boot$conf.int,
    beta = best.fit$beta,
    beta_CI = beta_CI$conf.interval,
    sigma_s = best.fit$sigma_s,
    mean_rup = best.fit$mean_rup,
    quality_summary = best.fit$quality_summary,
    fault_data = fault_df
  )
}


fault2michael <- function(x) {
  normals <- Fault_plane(x) |>
    Vec3() |>
    unclass()
  slips <- Fault_slip(x) |>
    Vec3() |>
    unclass()
  linear_stress_inversion(normals, slips)
}


.slip_inversion_michael <- function(x, friction = 0.6) {
  normals <- Fault_plane(x) |>
    Vec3() |>
    unclass()
  slips <- Fault_slip(x) |>
    Vec3() |>
    unclass()

  tau <- linear_stress_inversion(normals, slips)

  nx <- nrow(x)
  # tau0 <- tau / sqrt(sum(tau^2)) # normalize Frobenius norm

  tau_stress <- tau2stress(tau)
  sigma_vals <- tau_stress$sigma_vals
  principal_axes <- tau_stress$principal_axes


  # stress ratios:
  tau.params <- tau_params(tau)

  # Angles between the tangential traction predicted by the best stress tensor and the slip vector on each plane
  betas <- sapply(seq_len(nx), function(i) {
    int <- crossprod(Plane(principal_axes[2, ]), Plane(x[i, ])) |> Line()
    angle(int, Line(x[i, ]))
  }) #|> as.vector()
  betas <- ifelse(betas > 90, 180 - betas, betas)
  beta_mean <- tectonicr::circular_mean(betas, axial = FALSE)


  # Angle between slip planes and sigma 1
  theta <- vapply(seq_len(nx), function(i) {
    angle(Plane(x[i, ]), principal_axes[1, ])
  }, numeric(1))

  rup <- .rup(tau, normals, slips)
  quality <- ifelse(rup <= 50, "good",
                    ifelse(rup <= 75, "acceptable", "poor")
  )
  quality_summary <- c(n_good = sum(rup <= 50), n_acceptable = sum(rup > 50 & rup <= 75),
                       n_poor = sum(rup > 75))
  
  # Theoretically resolved shear stress on plane
  sigma_s_mean <- mean(shear_stress(sigma_vals[1], sigma_vals[3], theta))

  pr <- principal_fault(principal_axes[1, ], principal_axes[3, ], friction)

  shearnorm <- tau2shearnorm(tau, x, friction = friction)
  sigma_s <- shearnorm[, "shear"]
  sigma_n <- shearnorm[, "normal"]
  slip_tend <- slip_tendency(sigma_s, sigma_n)
  dilat_tend <- dilatation_tendency(sigma_vals[1], sigma_vals[3], sigma_n)


  SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(tau, symmetric = TRUE)$vectors),
    error = function(e) NA
  )

  list(
    stress_tensor = tau,
    principal_axes = principal_axes,
    principal_vals = sigma_vals,
    principal_faults = pr,
    SHmax = SHmax,
    R = tau.params$R,
    phi = tau.params$phi,
    bott = tau.params$bott,
    beta = beta_mean,
    mean_rup = mean(rup),
    quality_summary = quality_summary,
    sigma_s = sigma_s_mean,
    fault_data = data.frame(beta = betas, theta = theta, sigma_s = sigma_s, sigma_n = sigma_n, slip_tendency = slip_tend, dilational_tendency = dilat_tend, rup = rup, quality = quality)
  )
}


#' @importFrom MASS ginv
linear_stress_inversion <- function(n, s) {
  m <- nrow(n)

  # Eq 8:
  A <- do.call(
    rbind,
    lapply(seq_len(m), function(i) fault_normal_matrix(n[i, ]))
  )

  # Solve using generalized inverse (real part in case of rounding errors)

  # all produce same result:
  stress_vector <- Re(MASS::ginv(A) %*% as.vector(t(s)))
  # stress_vector <- qr.solve(A, as.vector(t(s)))
  # least_square_solv <- lm(as.vector(t(s)) ~ 0 + A)
  # stress_vector <- least_square_solv$coefficients

  # trace of stress tensor tau is usually assumed to be zero (Michael 1984):
  # trace(tau)  = sigma1 + sigma2 + sigma3 = 0 (Eq. 2)

  # Taken in combination with the constraint that the isotropic stress is zero, i.e. s33 = -(s11 + s22))
  stress_vector[6] <- -(stress_vector[1] + stress_vector[4])
  names(stress_vector) <- c("11", "12", "13", "22", "23", "33")

  stress_tensor <- matrix(c(
    stress_vector["11"], stress_vector["12"], stress_vector["13"],
    stress_vector["12"], stress_vector["22"], stress_vector["23"],
    stress_vector["13"], stress_vector["23"], stress_vector["33"]
  ), nrow = 3L, byrow = TRUE)

  return(stress_tensor)
}

fault_normal_matrix <- function(n) {
  # Eq. 8
  nsq <- n^2

  nmult <- -2 * n[1] * n[2] * n[3]

  A <- matrix(nrow = 5, ncol = 3)
  A[1, 1] <- n[1] * (nsq[2] + 2 * nsq[3])
  A[1, 2] <- n[2] * (-nsq[1] + nsq[3])
  A[1, 3] <- n[3] * (-2 * nsq[1] - nsq[2])

  A[2, 1] <- n[2] * (1 - 2 * nsq[1])
  A[2, 2] <- n[1] * (1 - 2 * nsq[2])
  A[2, 3] <- nmult

  A[3, 1] <- n[3] * (1 - 2 * nsq[1])
  A[3, 2] <- nmult
  A[3, 3] <- n[1] * (1 - 2 * nsq[3])

  A[4, 1] <- n[1] * (-nsq[2] + nsq[3])
  A[4, 2] <- n[2] * (nsq[1] + 2 * nsq[3])
  A[4, 3] <- n[3] * (-nsq[1] - 2 * nsq[2])

  A[5, 1] <- nmult
  A[5, 2] <- n[3] * (1 - 2 * nsq[2])
  A[5, 3] <- n[2] * (1 - 2 * nsq[3])

  t(A)
}



# Angelier (1990) --------------------------------------------------------------

#' Stress Inversion for Fault-Slip Data after Angelier (1990)
#'
#' Iterative direct inversion after the algorithm of Angelier (1990) and Mostafa (2005)
#'
#' @inheritParams slip_inversion_michael
#' @param max_iter integer. maximum Mostafa iteration count (default `50`)
#' @param tol numeric. convergence tolerance on max absolute change in TR elements between iterations (default `1e-6`)
#' @param n_psi integer. number of psi grid points for each step (default `361`)
#'
#' @details
#' The reduced stress tensor (Eq. 4.87) is parameterised as:
#' \deqn{T_R = \begin{bmatrix}
#'      \cos(\psi) & d & e \\
#'      d & \cos(\psi+2\pi/3) & f \\
#'      e & f & \cos(\psi + 4\pi/3)
#'    \end{bmatrix}
#'    }
#' with two normalisation constraints (Pascal, 2022; Eqs 4.88–4.89):
#'   1. \eqn{T_{11} + T_{22} + T_{33} = 0}            (deviator)
#'   2. \eqn{T_{11}^2 + T_{22}^2 + T_{33}^2 = 3/2}   (fixes \eqn{\lambda = \sqrt(3)/2})
#'
#' The four unknowns are \eqn{\psi}, \eqn{d}, \eqn{e}, \eqn{f}.
#'
#' Minimisation function (Eq. 4.101):
#'   \eqn{F_4 = \sum_i \upsilon_i^2},  where \eqn{\upsilon_i = \lambda * \hat{s}_i - \tau_i}
#'
#' \eqn{dF_4/d(d,e,f)} = 0 yields a 3x3 linear system in (\eqn{d},\eqn{e},\eqn{f}) given \eqn{\psi}.
#' \eqn{dF_4/d\psi = 0} is nonlinear; solved here by grid search + Brent refinement.
#'
#' Mostafa (2005) replaces the global \eqn{\lambda} with per-fault \eqn{\lambda_i} equal to
#' the shear traction magnitude on each plane and iterates until convergence.
#'
#' @references 
#' Angelier, J. (1990). Inversion of field data in fault tectonics to obtain the
#'  regional stress—III. A new rapid direct inversion method by analytical means. 
#'  Geophys. J. Int, 103, 363–376. <https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>
#'  
#' Pascal, C. (2022). Paleostress Inversion Techniques. Chapter 4, Sections 4.2.3 and 4.2.4.
#' 
#' Mostafa, M. E. (2005). Iterative direct inversion: An exact complementary 
#' solution for inverting fault-slip data to obtain palaeostresses. Computers 
#' & Geosciences, 31(8), 1059–1070. <https://doi.org/10.1016/j.cageo.2005.02.012>
#' 
#' @family stress-inversion
#'
#' @returns a named list
#' \describe{
#'   \item{`stress_tensor`}{matrix. Best-fit devitoric stress tensor in input coordinate frame}
#'   \item{`principal_axes`}{`"Line"` objects. Orientation of the principal stress axes as unit vectors (max to min)}
#'   \item{`tensor_params`}{the four tensor parameters (Eq. 4.87)}
#'   \item{`principal_vals`}{eigenvalues of the stress tensor (sigma1 >= sigma2 >= sigma3)}
#'   \item{`R`}{numeric. Stress shape ratio after Gephart & Forsyth (1984): \eqn{R = (\sigma_1 - \sigma_2)/(\sigma_1 - \sigma_3)}. Values ranging from 0 to 1, with 0 being
#'      \eqn{\sigma_1 = \sigma_2} and 1 being \eqn{\sigma_2 = \sigma_3}.}
#'  \item{`phi`}{numeric. Stress shape ratio after Angelier (1979): \eqn{\Phi = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)}. Values range between 0 (\eqn{\sigma_2 = \sigma_3}) and 1 (\eqn{\sigma_2 = \sigma_1}).}
#'  \item{`bott`}{numeric. Stress shape ratio after Bott (1959): \eqn{\R = (\sigma_3 - \sigma_1)/(\sigma_2 - \sigma_1)}. Values range between \eqn{-\infty} and \eqn{+\infty}.}
#'   \item{`alpha`}{numeric. The mean misfit angle in degrees. }
#'   \item{`mean_rup`}{mean RUP across all faults}
#'   \item{`quality_summary`}{counts per quality class}
#'  \item{`SHmax`}{numeric. Direction of maximum horizontal stress (in degrees)}
#'  \item{`beta`}{numeric. Average angle between the tangential traction predicted by the best stress tensor  and the slip vector on each plane. Should be close to 0.}
#'  \item{`sigma_s`}{numeric. Average resolved shear stress on each plane. Should be close to 1.}
#'  \item{`fault_data`}{`data.frame` containing the `alpha` angles, `beta` angles, i.e. the angles between sigma 1 and the plane normal,
#'  the resolved shear and normal stresses, the slip and dilation tendency on each plane, and the per-fault "ratio upsilon parameter" (RUP) estimator in percent, and the RUP-based quality.}
#'  \item{`n_iter`}{number of Mostafa iterations performed}
#'  }
#' @export
#'
#' @examples
#' # Use Angelier examples:
#' par(mfrow = c(1, length(angelier1990)))
#' 
#' # loop through dataset
#' invisible(lapply(angelier1990, function(x){
#' 
#' # inversion
#' res <- slip_inversion_angelier(x)
#' 
#' # stress shape
#' R_val <- round(res$R, 2)
#' 
#' # misfit
#' rup_val <- round(res$mean_rup, 2)
#'
#' # Plot the faults (color-coded by RUO%) and show the principal stress axes
#' stereoplot(title = "Iterative direct inversion", guides = FALSE)
#' stereo_shmax(res$SHmax)
#' fault_plot(x, col = assign_col(res$fault_data$rup))
#' points(res$principal_axes, col = 1:3, pch = 16, cex = 1.5)
#' text(res$principal_axes, label = rownames(res$principal_axes), 
#' col = 1:3, adj = -.25)
#' legend("topleft", col = 2:4, legend = rownames(res$principal_axes), pch = 16)
#' title(sub = bquote(R == .(R_val) ~ "|" ~ bar("RUP") == .(rup_val) * '%'))
#' }))
slip_inversion_angelier <- function(x,
                                    max_iter = 50L,
                                    tol = 1e-6,
                                    n_psi = 361L,
                                    friction = 0.6) {
  normals <- Fault_plane(x) |>
    Vec3() |>
    unclass()
  slips <- Fault_slip(x) |>
    Vec3() |>
    unclass()
  
  if (nrow(normals) < 4L) {
    stop("At least 4 fault slip measurements are required.")
  }
  
  N <- nrow(normals)
  lambda <- sqrt(3) / 2 # global lambda from normalised TR, Eq. 4.106
  
  # --- Step 1: initial Angelier (1990) inversion with global lambda ---
  res <- .angelier_step(normals, slips, lambda, n_psi)
  TR <- res$TR
  #
  # if (verbose)
  #   message(sprintf("Angelier init:  F4 = %.6f", res$F4))
  #
  # --- Step 2: Mostafa (2005) iteration ---
  # Replace global lambda with per-fault lambda_i = |tau_i| (shear traction
  # magnitude on each plane under the current TR) and re-invert.
  # Repeat until TR converges (max absolute element change < tol).
  
  for (iter in seq_len(max_iter)) {
    tau <- .shear_traction(TR, normals)
    tau_mag <- sqrt(rowSums(tau^2))
    
    # Fall back to global lambda for degenerate planes (zero shear)
    lambda_i <- ifelse(tau_mag > 1e-12, tau_mag, lambda)
    
    res_new <- .angelier_step(normals, slips, lambda_i, n_psi)
    TR_new <- res_new$TR
    
    delta <- max(abs(TR_new - TR))
    
    # if (verbose)
    #   message(sprintf("Mostafa iter %2d: F4 = %.6f  |dTR|_max = %.2e",
    #                   iter, res_new$F4, delta))
    
    TR <- TR_new
    res <- res_new
    
    if (delta < tol) break
  }
  
  # --- Step 3: Extract principal stresses ---
  # eig <- eigen(TR, symmetric = TRUE) # eigenvalues in DECREASING order
  # eigval <- eig$values
  # eigvec <- eig$vectors # columns are principal axes
  
  p <- tau2stress(TR)
  stress_shape <- tau_params(TR)
  
  # phi <- stress_shape$R
  
  # --- Step 4: Per-fault diagnostics ---
  tau_f <- .shear_traction(TR, normals)
  tau_norm <- sqrt(rowSums(tau_f^2))
  valid <- tau_norm > 1e-12
  
  tau_hat <- tau_f
  tau_hat[valid, ] <- tau_f[valid, ] / tau_norm[valid]
  
  # abs() handles the slip-sense ambiguity (Section 4.2.1, condition 4.16/4.17)
  cos_alpha <- pmax(-1, pmin(1, rowSums(tau_hat * slips)))
  alpha_deg <- acos(abs(cos_alpha)) * 180 / pi
  alpha_mean <- tectonicr::circular_mean(alpha_deg, axial = FALSE)
  
  rup <- .rup(TR, normals, slips, lambda)
  quality <- ifelse(rup <= 50, "good",
                    ifelse(rup <= 75, "acceptable", "poor")
  )
  quality_summary <- c(n_good = sum(rup <= 50), n_acceptable = sum(rup > 50 & rup <= 75),
                       n_poor = sum(rup > 75))
  
  
  tensor_params <- c(psi = res$psi, d = res$d, e = res$e, f = res$f)
  
  
  # Angles between the tangential traction predicted by the best stress tensor and the slip vector on each plane
  betas <- sapply(seq_len(N), function(i) {
    int <- crossprod(Plane(p$principal_axes[2, ]), Plane(x[i, ])) |> Line()
    angle(int, Line(x[i, ]))
  }) 
  betas <- ifelse(betas > 90, 180 - betas, betas)
  beta_mean <- tectonicr::circular_mean(betas, axial = FALSE)
  
  
  # Angle between slip planes and sigma 1
  theta <- vapply(seq_len(N), function(i) {
    angle(Plane(x[i, ]), p$principal_axes[1, ])
  }, numeric(1))
  
  # Theoretically resolved shear stress on plane
  sigma_s_mean <- mean(shear_stress(p$sigma_vals[1], p$sigma_vals[3], theta))
  
  shearnorm <- tau2shearnorm(TR, x, friction = friction)
  sigma_s <- shearnorm[, "shear"]
  sigma_n <- shearnorm[, "normal"]
  slip_tend <- slip_tendency(sigma_s, sigma_n)
  dilat_tend <- dilatation_tendency(p$sigma_vals[1], p$sigma_vals[3], sigma_n)
  
  SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(TR, symmetric = TRUE)$vectors),
    error = function(e) {
      SH(p$principal_axes[1, ], p$principal_axes[2, ], p$principal_axes[3, ], R = stress_shape$R)
    }
  )
  
  pfaults <- principal_fault(p$principal_axes[1, ], p$principal_axes[3, ], friction)
  
  list(
    stress_tensor = TR,
    tensor_params = tensor_params,
    principal_axes = p$principal_axes,
    principal_vals = p$sigma_vals,
    principal_faults = pfaults,
    R = stress_shape$R,
    phi = stress_shape$phi,
    bott = stress_shape$bott,
    sigma_s = sigma_s_mean,
    alpha = alpha_mean,
    beta = beta_mean,
    mean_rup = mean(rup),
    quality_summary = quality_summary,
    n_iter = iter,
    SHmax = SHmax,
    fault_data = data.frame(alpha = alpha_deg, beta = betas, theta = theta, sigma_s = sigma_s, sigma_n = sigma_n, slip_tendency = slip_tend, dilational_tendency = dilat_tend, rup = rup, quality = quality)
  )
}


# Build the 3x3 reduced stress tensor from the four unknowns (Eq. 4.87)
.build_TR <- function(psi, d, e, f) {
  cp <- cos(psi)
  sp <- sin(psi)
  r3 <- sqrt(3) / 2
  
  # Diagonal: T11=cos(psi), T22=cos(psi+2pi/3), T33=cos(psi+4pi/3)
  # Using addition formula: cos(psi+2pi/3) = -cp/2 - sp*r3
  #                         cos(psi+4pi/3) = -cp/2 + sp*r3
  matrix(
    c(
      cp, d, e,
      d, -cp / 2 - sp * r3, f,
      e, f, -cp / 2 + sp * r3
    ),
    nrow = 3L, byrow = FALSE
  ) # column-major: rows of byrow=F are columns
}

# Shear traction vectors for N fault planes (returns N x 3 matrix)
.shear_traction <- function(TR, normals) {
  tractions <- normals %*% TR # N x 3
  normal_proj <- rowSums(tractions * normals) # N (scalar normal stress)
  tractions - normal_proj * normals # N x 3 shear vectors
}

# RUP estimator per fault, Eq. 4.107:
#   RUP_i = 100 * |upsilon_i| / (sqrt(3)/2)
# where upsilon_i = lambda * s_hat_i - tau_i  (Eq. 4.100)
.rup <- function(TR, normals, slips, lambda = sqrt(3) / 2) {
  tau <- .shear_traction(TR, normals)
  ups_m <- sqrt(rowSums((slips * lambda - tau)^2))
  100 * ups_m / (sqrt(3) / 2)
}


# One Angelier (1990) direct inversion step
#
# For a given lambda vector (scalar or length-N for Mostafa mode) this
# function:
#   1. Accumulates the coefficient matrices for the 3x3 linear system in
#      (d, e, f) as a function of psi (from dF4/d{d,e,f} = 0).
#   2. Evaluates F4 at a grid of psi values by solving the 3x3 system
#      at each node.
#   3. Refines the minimum with Brent's method.
#   4. Returns the optimal (psi, d, e, f) and the corresponding TR.
#' @importFrom stats optimise
.angelier_step <- function(normals, slips, lambdas, n_psi = 181L) {
  N <- nrow(normals)
  lv <- if (length(lambdas) == 1L) rep(lambdas, N) else lambdas
  
  n1 <- normals[, 1L]
  n2 <- normals[, 2L]
  n3 <- normals[, 3L]
  s1 <- slips[, 1L]
  s2 <- slips[, 2L]
  s3 <- slips[, 3L]
  r3 <- sqrt(3) / 2
  
  # Symbolic expansion of tau_i components in terms of (cp, sp, d, e, f)
  #
  # TR diagonal in terms of cp = cos(psi), sp = sin(psi):
  #   T11 =  cp
  #   T22 = -cp/2 - sp*r3
  #   T33 = -cp/2 + sp*r3
  #
  # Stress vector on plane i:
  #   sigma_ix = cp*n1       + d*n2  + e*n3
  #   sigma_iy = d*n1        + (-cp/2 - sp*r3)*n2  + f*n3
  #   sigma_iz = e*n1        + f*n2  + (-cp/2 + sp*r3)*n3
  #
  # Normal stress (scalar):
  #   sigma_n = cp*Ac + sp*As + 2*d*n1*n2 + 2*e*n1*n3 + 2*f*n2*n3
  # where:
  #   Ac = n1^2 - n2^2/2 - n3^2/2
  #   As = r3*(n3^2 - n2^2)
  #
  # Shear vector: tau_i = sigma_i - sigma_n * n_i
  # Components (linear in cp, sp, d, e, f):
  #
  #   tau_ix = cp*(n1 - Ac*n1)   + sp*(-As*n1)
  #           + d*(n2 - 2*n1^2*n2) + e*(n3 - 2*n1^2*n3) + f*(-2*n1*n2*n3)
  #
  #   tau_iy = cp*(-n2/2 - Ac*n2) + sp*(-r3*n2 - As*n2)
  #           + d*(n1 - 2*n1*n2^2) + e*(-2*n1*n2*n3)    + f*(n3 - 2*n2^2*n3)
  #
  #   tau_iz = cp*(-n3/2 - Ac*n3) + sp*(r3*n3 - As*n3)
  #           + d*(-2*n1*n2*n3)   + e*(n1 - 2*n1*n3^2)  + f*(n2 - 2*n2*n3^2)
  
  Ac <- n1^2 - n2^2 / 2 - n3^2 / 2
  As <- r3 * (n3^2 - n2^2)
  
  # Coefficient arrays: tau_ix = Px_c*cp + Px_s*sp + Px_d*d + Px_e*e + Px_f*f
  Px_c <- n1 * (1 - Ac)
  Px_s <- -As * n1
  Px_d <- n2 * (1 - 2 * n1^2)
  Px_e <- n3 * (1 - 2 * n1^2)
  Px_f <- -2 * n1 * n2 * n3
  
  Py_c <- -n2 * (0.5 + Ac)
  Py_s <- -n2 * (r3 + As)
  Py_d <- n1 * (1 - 2 * n2^2)
  Py_e <- -2 * n1 * n2 * n3
  Py_f <- n3 * (1 - 2 * n2^2)
  
  Pz_c <- -n3 * (0.5 + Ac)
  Pz_s <- n3 * (r3 - As)
  Pz_d <- -2 * n1 * n2 * n3
  Pz_e <- n1 * (1 - 2 * n3^2)
  Pz_f <- n2 * (1 - 2 * n3^2)
  
  # dF4/dd = 0:  sum[ d/dd(tau^2) - 2*li * d/dd(s.sigma) ] = 0
  #
  # d/dd(tau_i^2) = 2 * (dTR/dd * n_i) . tau_i
  #   dTR/dd * n_i = (n2, n1, 0)^T
  #   => 2 * (n2*tau_ix + n1*tau_iy)
  #
  # d/dd(s_i.sigma_i) = s_i . (dTR/dd * n_i) = s1*n2 + s2*n1
  #
  # Equation for d:  A*d + D*e + E*f = Gc*cp + Gs*sp + U
  #   A  = sum(n2*Px_d + n1*Py_d)
  #   D  = sum(n2*Px_e + n1*Py_e)   [= sum(n3*Px_d + n1*Pz_d) by symmetry]
  #   E  = sum(n2*Px_f + n1*Py_f)
  #   Gc = sum(n2*Px_c + n1*Py_c)
  #   Gs = sum(n2*Px_s + n1*Py_s)
  #   U  = sum(lv*(s1*n2 + s2*n1))
  #
  # Similarly for dF4/de (lever = (n3,0,n1)) and dF4/df (lever = (0,n3,n2)).
  
  A_ <- sum(n2 * Px_d + n1 * Py_d)
  D_ <- sum(n2 * Px_e + n1 * Py_e)
  E_ <- sum(n2 * Px_f + n1 * Py_f)
  Gc_ <- sum(n2 * Px_c + n1 * Py_c)
  Gs_ <- sum(n2 * Px_s + n1 * Py_s)
  U_ <- sum(lv * (s1 * n2 + s2 * n1))
  
  B_ <- sum(n3 * Px_e + n1 * Pz_e)
  F_ <- sum(n3 * Px_f + n1 * Pz_f)
  Hc_ <- sum(n3 * Px_c + n1 * Pz_c)
  Hs_ <- sum(n3 * Px_s + n1 * Pz_s)
  V_ <- sum(lv * (s1 * n3 + s3 * n1))
  
  C_ <- sum(n3 * Py_f + n2 * Pz_f)
  Ic_ <- sum(n3 * Py_c + n2 * Pz_c)
  Is_ <- sum(n3 * Py_s + n2 * Pz_s)
  W_ <- sum(lv * (s2 * n3 + s3 * n2))
  
  # Symmetric 3x3 coefficient matrix (constant across psi)
  lhs_mat <- matrix(c(
    A_, D_, E_,
    D_, B_, F_,
    E_, F_, C_
  ), 3L, 3L)
  
  # Solve for (d, e, f) given psi
  .solve_def <- function(psi_val) {
    cp <- cos(psi_val)
    sp <- sin(psi_val)
    rhs <- c(
      U_ - Gc_ * cp - Gs_ * sp,
      V_ - Hc_ * cp - Hs_ * sp,
      W_ - Ic_ * cp - Is_ * sp
    )
    tryCatch(solve(lhs_mat, rhs), error = function(e) rep(NA_real_, 3L))
  }
  
  # F4 evaluated at a given psi (after solving for d,e,f)
  .F4 <- function(psi_val) {
    def <- .solve_def(psi_val)
    if (anyNA(def)) {
      return(.Machine$double.xmax)
    }
    TR <- .build_TR(psi_val, def[1L], def[2L], def[3L])
    tau <- .shear_traction(TR, normals) # N x 3
    # F4 = sum[ lv^2 + |tau|^2 - 2*lv*(s.tau) ]  (Eq. 4.103)
    sum(lv^2 + rowSums(tau^2) - 2 * lv * rowSums(slips * tau))
  }
  
  # Grid search over psi in [0, pi)  (F4 has period pi in psi)
  psi_grid <- seq(0, pi, length.out = n_psi)
  F4_grid <- vapply(psi_grid, .F4, numeric(1L))
  best_idx <- which.min(F4_grid)
  best_psi <- psi_grid[best_idx]
  
  # Brent refinement within one grid step of the minimum
  step <- pi / (n_psi - 1L)
  opt <- tryCatch(
    stats::optimise(.F4, interval = c(
      max(0, best_psi - step),
      min(pi, best_psi + step)
    )),
    error = function(e) list(minimum = best_psi, objective = F4_grid[best_idx])
  )
  
  psi_opt <- opt$minimum
  def_opt <- .solve_def(psi_opt)
  TR_opt <- .build_TR(psi_opt, def_opt[1L], def_opt[2L], def_opt[3L])
  
  list(
    psi = psi_opt,
    d = def_opt[1L],
    e = def_opt[2L],
    f = def_opt[3L],
    TR = TR_opt,
    F4 = opt$objective
  )
}


# PT technique ------------------------------------------------------------------


#' Simple fault analysis
#'
#' The PT-techniques is a graphical solution of the *Wallace-Bott hypothesis*, i.e. fault slip occurs parallel to the maximum shear stress.
#' It calculates PT-axes, kinematic planes (also movement planes), and the dihedra separation plane.
#'
#' @inheritParams slip_inversion
#' @param ptangle numeric. angle between P and T axes in degrees (90&deg; by default).
#'
#' @returns list. `p` and `t` are the P and T axes as `"Line"` objects,
#' `m` and `d` are the M-planes and the dihedra separation planes as `"Plane"` objects
#' @export
#'
#' @family stress-inversion
#'
#' @examples
#' par(mfrow = c(1, length(angelier1990)))
#' invisible(lapply(angelier1990, function(x){ 
#'   xpt <- Fault_PT(x)
#'   
#'   stereoplot(guides = FALSE)
#'   angelier(x, col = 'grey')
#'   points(xpt$p, pch = 16, cex = 1.5, col = 1)
#'   points(xpt$t, pch = 16, cex = 1.5, col = 2)
#'   stereo_confidence(xpt$p, pch = 16, cex = 1.5, col = 1)
#'   stereo_confidence(xpt$t, pch = 16, cex = 1.5, col = 2)
#'   }))
Fault_PT <- function(x, ptangle = 90) {
  ptangle <- rep(deg2rad(ptangle), nrow(x))
  x_corr <- misfit_pair(x)
  xp <- x_corr$fvec
  xl <- x_corr$lvec * ifelse(x[, 5] < 0, -1, 1)

  pvec <- rotate(xp, crossprod(xl, xp), -ptangle / 2)
  tvec <- rotate(xp, crossprod(xl, xp), ptangle / 2)

  p <- Line(pvec)
  t <- Line(tvec)

  mvec <- crossprod(xl, xp)
  m <- Plane(mvec)

  dvec <- crossprod(crossprod(xl, xp), xp)
  d <- Plane(dvec)

  list(p = p, t = t, m = m, d = d)
}

rot_mean <- function(x){
  x |> pair2rot() |> rotMean() |> as.Rotation() |> rot2pair(fault = TRUE)
}

#' Simple Statistical Fault-Slip Inversion
#' 
#' Measurements of fault-slip data are often scattered due to measurement 
#' errors and the wavy nature of fault planes and fault striations/slickenlines. 
#' The fault scatter is large due to noise, rather than representing the actual 
#' geometry of
#' the fault set. The idea of this algorithm is to cluster the fault data set to 
#' identify the conjugate set of 
#' faults and their the mean orientation Using the Wallace-Bott Hypothesis 
#' and Anderson's theory, it then calculates the orientation of the principal stresses, 
#' and uses the angles to the fault planes to derive the a best-fit stress shape parameter R. 
#'
#' @param x object of class `"Fault"`
#' @param cluster_fun function for cluster, must have number of desired cluster 
#' as second input outputs and vector `cluster`. The default is [stats::kmeans()]
#' @param n_grid integer. Number to optimize grid search for stress shape parameter R
#' 
#' @family stress-inversion
#'
#' @returns list.
#' @export
#'
#' @examples
#' tym <- slip_inversion_simple(angelier1990$TYM)
#' stereoplot(title = 'TYM', sub = paste0('beta: ', round(tym$beta, 2), 
#' " deg | R: ", round(tym$R, 2)))
#' hoeppener(angelier1990$TYM, col = assign_col(tym$beta_angles))
#' angelier(tym$mean_planes, pch = 16, col = viridis::magma(2, end = 0.8), cex = 1)
#' points(tym$principal_axes, pch = 16, col = viridis::rocket(3, end = 0.8), cex = 1)
#' text(tym$principal_axes, labels = rownames(tym$principal_axes), 
#' col = viridis::rocket(3, end = 0.8), cex = 1, adj = c(-.25, -.25))
#' 
#' avb <- slip_inversion_simple(angelier1990$AVB)
#' stereoplot(title = 'AVB', sub = paste0('beta: ', round(avb$beta, 2), 
#' " deg | R: ", round(avb$R, 2)))
#' hoeppener(angelier1990$AVB, col = assign_col(avb$beta_angles))
#' angelier(avb$mean_planes, pch = 16, col = viridis::magma(2, end = 0.8), cex = 1)
#' points(avb$principal_axes, pch = 16, col = viridis::rocket(3, end = 0.8), cex = 1)
#' text(avb$principal_axes, labels = rownames(avb$principal_axes), 
#' col = viridis::rocket(3, end = 0.8), cex = 1, adj = c(-.25, -.25))
slip_inversion_simple <- function(x, cluster_fun = stats::kmeans, n_grid = 1000L){
  td <- dist.Pair(x)
  tdcluster <- cluster_fun(td, 2)$cluster
  
  tsplit <- split.data.frame(unclass(x), tdcluster) |> 
    lapply(as.Fault)
  
  tclustermeans <- lapply(tsplit, rot_mean)
  tclustermeansP <- lapply(tclustermeans, Plane)
  
  # The intermediate stress axis (σ₂) represents the direction of no displacement within the fault plane.
  sigma2 <- crossprod(tclustermeansP[[1]], tclustermeansP[[2]])
  
  r <- which.max(table(tsplit[[1]][, 'sense']))
  ang <- angle(tclustermeansP[[1]], tclustermeansP[[2]])/2
  sigma1 <- rotate(tclustermeansP[[1]], sigma2, r*ang)
  
  if(!.near(angle(sigma1, tclustermeansP[[1]]), angle(sigma1, tclustermeansP[[2]]))) sigma1 <- rotate(tclustermeansP[[1]], sigma2, -r*ang)
  
  sigma3 <- crossprod(sigma1, sigma2)
  
  principal_axes <- rbind(sigma1, sigma2, sigma3)
  rownames(principal_axes) <- c('sigma1', 'sigma2', 'sigma3')
  
  # Angles between the tangential traction predicted by the best stress tensor and the slip vector on each plane
  betas <- sapply(seq_len(nrow(x)), function(i) {
    int <- crossprod(Plane(principal_axes[2, ]), Plane(x[i, ])) |> Line()
    angle(int, Line(x[i, ]))
  }) 
  betas <- ifelse(betas > 90, 180 - betas, betas)
  beta_mean <- tectonicr::circular_mean(betas, axial = FALSE)
  
  # Angle between slip planes and sigma 1
  theta <- sapply(seq_len(nrow(x)), function(i) {
    angle(Plane(x[i, ]), principal_axes[1, ])
  })
  
  res <- list(principal_axes = principal_axes, beta = beta_mean, theta = theta, x_split = tsplit, mean_planes = do.call(rbind, tclustermeans), cluster = tdcluster, beta_angles = betas)
  #Rlist <- .solve_R(stress_eigenvecs = principal_axes, faults = x, n_grid)
  Rlist <- list(R = NA)
  append(res, Rlist)
}


.solve_R <- function(stress_eigenvecs, faults, n_grid = 1000L) {
  stress_eigenvecs = unclass(Vec3(stress_eigenvecs))
  fault_normals = unclass(Vec3(Plane(faults)))
  slip_dirs = unclass(Vec3(Ray(faults)))
  
  e1 <- stress_eigenvecs[1L, ]
  e2 <- stress_eigenvecs[2L, ]
  
  misfit <- function(R) {
    sigma <- tcrossprod(e1) + R * tcrossprod(e2)
    tractions    <- fault_normals %*% t(sigma)
    normals_proj <- rowSums(tractions * fault_normals)
    shear        <- tractions - normals_proj * fault_normals
    shear_hat    <- shear / sqrt(rowSums(shear^2))
    cos_alpha    <- pmax(-1, pmin(1, rowSums(shear_hat * slip_dirs)))
    mean(acos(cos_alpha)^2)
  }
  
  R_grid   <- seq(0, 1, length.out = n_grid)
  misfits  <- vapply(R_grid, misfit, numeric(1L))
  best_idx <- which.min(misfits)
  step     <- 1 / (n_grid - 1L)
  lower    <- max(0, R_grid[best_idx] - step)
  upper    <- min(1, R_grid[best_idx] + step)
  opt      <- stats::optimise(misfit, interval = c(lower, upper))
  
  list(
    R         = opt$minimum,
    misfit_deg = sqrt(opt$objective) * 180 / pi,
    R_grid    = R_grid[best_idx],
    converged = abs(opt$minimum - R_grid[best_idx]) < step
  )
}

# Stress tensor ----------------------------------------------------------------

fault_instability_criterion <- function(x, R, mu = 0.6) {
  n <- Plane(x) |> Vec3()
  # n: normal to plane
  # R: stress shape ratio
  # mu: friction (e.g. Byerlee: 0.6 - 0.85)

  # Instability ranges from 0 (most stable) to 1 (most unstable)
  # The most unstable fault is the optimally oriented fault for shear faulting

  # sigma1 <- 1
  # sigma2 <- 1 - 2 * R
  # sigm3 <- -1

  nsq <- n^2
  sigma <- nsq[, 1] + (1 - 2 * R) * nsq[, 2] - nsq[, 3]
  tau <- nsq[, 1] + (1 - R)^2 * nsq[, 2] + nsq[, 3] - (nsq[, 1] + (1 - 2 * R) * nsq[, 2] - nsq[, 3])

  tau_c <- 1 / (sqrt(1 + mu^2))
  sigma_c <- -mu / (1 + mu^2)

  instability <- (tau - mu * (sigma - 1)) / (mu + sqrt(1 + mu^2))
  return(instability)
}


#' Calculate Principal Fault planes from stress vectors and friction
#'
#' @param s1,s3 Principal stress vectors as spherical objects
#' @param friction numeric. Coefficient of friction
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' res_TYM <- slip_inversion(angelier1990$TYM, n_iter = 10)
#' pr_TYM <- principal_fault(res_TYM$principal_axes[1, ], res_TYM$principal_axes[3, ])
#'
#' stereoplot()
#' fault_plot(angelier1990$TYM, col = "grey")
#' fault_plot(pr_TYM, col = "red")
#' points(res_TYM$principal_axes, pch = 16)
principal_fault <- function(s1, s3, friction = 0.6) {
  mu <- 0.5 * atan(1 / friction)

  s1_vec <- Vec3(s1)
  s3_vec <- Vec3(s3)

  n1 <- sin(mu) * s1_vec - cos(mu) * s3_vec
  u1 <- cos(mu) * s1_vec + sin(mu) * s3_vec
  sense1 <- ifelse(u1[, 3] >= 0, 1, -1)

  n2 <- sin(-mu) * s1_vec - cos(-mu) * s3_vec
  u2 <- cos(-mu) * s1_vec + sin(-mu) * s3_vec
  sense2 <- ifelse(u2[, 3] >= 0, 1, -1)

  n <- rbind(n1, n2)
  u <- rbind(u1, u2)
  sense <- c(sense1, sense2)

  Fault(Plane(n), Line(u), sense = sense)
}

# Extract principal stress from stress tensor
tau2stress <- function(tau) {
  eig <- eigen(tau, symmetric = TRUE)
  sigma_vals <- eig$values

  principal_axes <- t(eig$vectors) |>
    as.Vec3() |>
    Line() # sigma1, sigma2, sigma3
  names(sigma_vals) <- rownames(principal_axes) <- c("sigma1", "sigma2", "sigma3")

  list(sigma_vals = sigma_vals, principal_axes = principal_axes)
}

# Calculate normal and shear stress components for faults
tau2shearnorm <- function(tau, fault, friction) {
  # Principal fault planes
  stess <- tau2stress(tau)
  np <- principal_fault(s1 = stess$principal_axes[1, ], s3 = stess$principal_axes[3, ], friction = friction)
  np1 <- Plane(np[1, ]) |> Vec3()
  np2 <- Plane(np[2, ]) |> Vec3()

  # Fault normals in Cartesian coordinates
  n <- Plane(fault) |> Vec3()
  n1 <- n[, 1]
  n2 <- n[, 2]
  n3 <- n[, 3]

  # Shear and normal stress
  tau_normal <- tau[1, 1] * n1^2 + tau[2, 2] * n2^2 + tau[3, 3] * n3^2 +
    2 * (tau[1, 2] * n1 * n2 + tau[1, 3] * n1 * n3 + tau[2, 3] * n2 * n3)

  total <- tau %*% rbind(n1, n2, n3)
  tau_total <- sqrt(colSums(total^2))
  tau_shear <- sqrt(tau_total^2 - tau_normal^2)

  # Fault–principal deviation and half-plane
  # dot1 <- abs(n1 * np1[1] + n2 * np2[1] + n3 * np3[1])
  # dot2 <- abs(n1 * np1[2] + n2 * np2[2] + n3 * np3[2])
  #
  # dev1 <- acos(pmin(dot1, 1)) * 180 / pi
  # dev2 <- acos(pmin(dot2, 1)) * 180 / pi
  dev1 <- angle(n, np1)
  dev2 <- angle(n, np2)

  tau_shear <- ifelse(dev1 < dev2, tau_shear, -tau_shear)

  cbind(normal = tau_normal, shear = tau_shear)
}

tau_params <- function(tau) {
  tau_stress <- tau2stress(tau)
  sigma_vals <- tau_stress$sigma_vals
  principal_axes <- tau_stress$principal_axess

  # stress ratios:
  R <- (sigma_vals[1] - sigma_vals[2]) / (sigma_vals[1] - sigma_vals[3]) # Gephart & Forsyth 1984
  phi <- (sigma_vals[2] - sigma_vals[3]) / (sigma_vals[1] - sigma_vals[3]) # Angelier 1979
  shape_ratio_bott <- (sigma_vals[3] - sigma_vals[1]) / (sigma_vals[2] - sigma_vals[1]) # Bott, Simon-Gomez
  list(R = R, phi = phi, bott = shape_ratio_bott)
}