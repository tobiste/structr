#' Stress Inversion for Fault-Slip Data
#'
#' Linear stress inversion (based on Michael, 1984) determines the orientation
#' of the principal stresses from fault slip data.
#' Confidence intervals are estimated by bootstrapping.
#' This inversion is simplified by the assumption that the magnitude of the
#' tangential traction on the various fault planes, at the time of rupture, is similar.
#'
#' @param x `"Fault"` object where the rows are the observations, and the columns the coordinates.
#' @param boot integer. Number of bootstrap samples (10 by default)
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
#'
#' @importFrom stats t.test
#'
#' @seealso [Fault_PT()] for a simple P-T stress analysis,
#'  [SH()] and [SH_from_tensor()] to calculate the azimuth of the maximum horizontal stress;
#'  [Mohr_plot()] for graphical representation of the deviatoric stress tensor.
#'
#' @examples
#' set.seed(20250411)
#' # Use Angelier examples:
#' res_TYM <- slip_inversion(angelier1990$TYM, boot = 100, n = 1000, res = 100)
#'
#' # Plot the faults (color-coded by beta angle) and show the principal stress axes
#' stereoplot(title = "Tymbaki, Crete, Greece", guides = FALSE)
#' fault_plot(angelier1990$TYM, col = "gray80")
#' stereo_confidence(res_TYM$principal_axes_conf$sigma1, col = 2)
#' stereo_confidence(res_TYM$principal_axes_conf$sigma2, col = 3)
#' stereo_confidence(res_TYM$principal_axes_conf$sigma3, col = 4)
#' text(res_TYM$principal_axes, label = rownames(res_TYM$principal_axes), col = 2:4, adj = -.25)
#' legend("topleft", col = 2:4, legend = rownames(res_TYM$principal_axes), pch = 16)
#'
#' res_AVB <- slip_inversion(angelier1990$AVB, boot = 100, n = 1000, res = 100)
#' stereoplot(title = "Agia Varvara, Crete, Greece", guides = FALSE)
#' fault_plot(angelier1990$AVB, col = "gray80")
#' stereo_confidence(res_AVB$principal_axes_conf$sigma1, col = 2)
#' stereo_confidence(res_AVB$principal_axes_conf$sigma2, col = 3)
#' stereo_confidence(res_AVB$principal_axes_conf$sigma3, col = 4)
#' text(res_AVB$principal_axes, label = rownames(res_AVB$principal_axes), col = 2:4, adj = -.25)
#' legend("topleft", col = 2:4, legend = rownames(res_AVB$principal_axes), pch = 16)
slip_inversion <- function(x, boot = 100L, conf.level = 0.95, friction = 0.6, ...) {
  best.fit <- slip_inversion0(x, friction)
  fault_df <- best.fit$fault_data
  nx <- nrow(x)

  # bootstrap results
  boot_results <- lapply(1:boot, function(i) {
    idx <- sample.int(nx, replace = TRUE)
    x_sample <- x[idx, ]
    slip_inversion0(x_sample)
  })

  # calculate confidence intervals from bootstrap results
  sigma_vec1 <- do.call(rbind, lapply(boot_results, function(x) {
    x$principal_axes[1, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)

  sigma_vec2 <- do.call(rbind, lapply(boot_results, function(x) {
    x$principal_axes[2, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)

  sigma_vec3 <- do.call(rbind, lapply(boot_results, function(x) {
    x$principal_axes[3, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, ...)

  R_boot <- vapply(boot_results, function(x) {
    x$R
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  phi_boot <- vapply(boot_results, function(x) {
    x$phi
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  bott_boot <- vapply(boot_results, function(x) {
    x$bott
  }, FUN.VALUE = numeric(1)) |>
    stats::t.test(conf.level = conf.level)

  sigma_boot0 <- vapply(boot_results, function(x) {
    x$principal_vals
  }, FUN.VALUE = numeric(3)) |> t()
  sigma_boot <- sapply(1:3, function(col) {
    sigma_boot_col <- stats::t.test(sigma_boot0[, col], conf.level = conf.level)
    sigma_boot_col$conf.int
  })
  colnames(sigma_boot) <- names(best.fit$principal_vals)

  beta_CI <- tectonicr::confidence_interval(fault_df$beta, conf.level = conf.level, axial = FALSE)

  SHmax_CI <- vapply(boot_results, function(x) {
    x$SHmax
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
    fault_data = fault_df
  )
}

slip_inversion0 <- function(x, friction = 0.6) {
  tau <- linear_stress_inversion(x)
  # tau0 <- tau / sqrt(sum(tau^2)) # normalize Frobenius norm

  tau_stress <- tau2stress(tau)
  sigma_vals <- tau_stress$sigma_vals
  principal_axes <- tau_stress$principal_axes


  # stress ratios:
  R <- (sigma_vals[1] - sigma_vals[2]) / (sigma_vals[1] - sigma_vals[3]) # Gephart & Forsyth 1984
  phi <- (sigma_vals[2] - sigma_vals[3]) / (sigma_vals[1] - sigma_vals[3]) # Angelier 1979
  shape_ratio_bott <- (sigma_vals[3] - sigma_vals[1]) / (sigma_vals[2] - sigma_vals[1]) # Bott, Simon-Gomez

  # Angles between the tangential traction predicted by the best stress tensor and the slip vector on each plane
  betas <- sapply(seq_len(nrow(x)), function(i) {
    int <- crossprod(Plane(principal_axes[2, ]), Plane(x[i, ])) |> Line()
    angle(int, Line(x[i, ]))
  }) #|> as.vector()
  betas <- ifelse(betas > 90, 180 - betas, betas)
  beta_mean <- tectonicr::circular_mean(betas, axial = FALSE)


  # Angle between slip planes and sigma 1
  theta <- sapply(seq_len(nrow(x)), function(i) {
    angle(Plane(x[i, ]), principal_axes[1, ])
  })

  # Theoretically resolved shear stress on plane
  sigma_s_mean <- mean(shear_stress(sigma_vals[1], sigma_vals[3], theta))

  pr <- principal_fault(principal_axes[1, ], principal_axes[3, ], friction)

  shearnorm <- tau2shearnorm(tau, x, friction = friction)
  sigma_s <- shearnorm[, "shear"]
  sigma_n <- shearnorm[, "normal"]
  slip_tend <- slip_tendency(sigma_s, sigma_n)
  dilat_tend <- dilatation_tendency(sigma_vals[1], sigma_vals[3], sigma_n)


  SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(tau)$vectors),
    error = function(e) NA
  )

  list(
    stress_tensor = tau,
    principal_axes = principal_axes,
    principal_vals = sigma_vals,
    principal_faults = pr,
    SHmax = SHmax,
    R = unname(R),
    phi = unname(phi),
    bott = unname(shape_ratio_bott),
    beta = beta_mean,
    sigma_s = sigma_s_mean,
    fault_data = data.frame(beta = betas, theta = theta, sigma_s = sigma_s, sigma_n = sigma_n, slip_tendency = slip_tend, dilational_tendency = dilat_tend)
  )
}


#' @importFrom MASS ginv
linear_stress_inversion <- function(x) {
  m <- nrow(x)

  # plane normal
  n <- Plane(x) |>
    Vec3() |>
    unclass()

  # slip vector
  s <- Ray(x) |>
    Vec3() |>
    unclass()
  # rake <- structr::Fault_rake(fault) * pi / 180

  # Eq 8:
  A <- do.call(
    rbind,
    lapply(1:m, function(i) fault_normal_matrix(n[i, ]))
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
  ), nrow = 3, byrow = TRUE)

  return(stress_tensor)
}

fault_normal_matrix <- function(n) {
  # Eq. 8
  nsq <- n^2

  nmult <- -2 * n[1] * n[2] * n[3]

  # A <- matrix(nrow = 3, ncol = 5)
  # A[1, 1] = n[1] - n[1]^3 + n[1] * n[3]^2
  # A[1, 2] = n[2] - 2 * n[2] * n[1]^2
  # A[1, 3] = n[3] - 2 * n[3] * n[1]^2
  # A[1, 4] = -n[1] * n[2]^2 + n[1] * n[3]^2
  # A[1, 5] = nmult
  #
  # A[2, 1] = -n[2] * n[1]^2 + n[2] * n[3]^2
  # A[2, 2] = n[1] - 2 * n[1] * n[2]^2
  # A[2, 3] = nmult
  # A[2, 4] = n[2] - 2 * n[2]^3 + n[2] * n[3]^2
  # A[2, 5] = n[3] - 2 * n[3]* n[2]^2
  #
  # A[3, 1] = -n[3] * n[1]^2 + n[3] + n[3]^2
  # A[3, 2] = nmult
  # A[3, 3] = n[1] - 2 * n[1] * n[3]^2
  # A[3, 4] = -n[2]^2 * n[3] - n[3] + n[3]^3
  # A[3, 5] = n[2] - 2*n[2] * n[3]^2
  #
  # A


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
#' res_TYM <- slip_inversion(angelier1990$TYM, boot = 10)
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
  eig <- eigen(tau)
  # sigma_vals <- sort(eig$values, decreasing  = TRUE)
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
#' @seealso [slip_inversion()]
#'
#' @examples
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' Fault_PT(f)
Fault_PT <- function(x, ptangle = 90) {
  ptangle <- deg2rad(ptangle)
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

#' Orthogonalization of plane and line measurement
#'
#' Both the line and the plane  are rotated in opposite directions by half the angle between
#' the slip and the plane normal vector about the vector that is perpendicular to both.
#'
#' @param x object of class `"Pair"` or `"Fault"`
#'
#' @returns `misfit_pair` returns a list with the orthogonalized plane and line measurements
#' (as 3d vectors) and the misfit angles (in radians). `correct_pair` returns the orthogonalized vectors.
#' @name pair_correct
#' @examples
#' p <- Pair(120, 60, 110, 58, correction = FALSE)
#' misfit_pair(p)
#'
#' correct_pair(p)
NULL

#' @rdname pair_correct
#' @export
misfit_pair <- function(x) {
  stopifnot(is.Pair(x))
  fvec <- Plane(x[, 1], x[, 2]) |> Vec3()
  lvec <- Line(x[, 3], x[, 4]) |> Vec3()

  misfit <- abs(pi / 2 - angle(fvec, lvec))
  # for (i in 1:length(misfit)) {
  #   if (misfit[i] > deg2rad(20)) {
  #     warning(paste("Mifit angle is", misfit[i], "degrees."))
  #   }
  # }

  # Warn if misfit > 20
  exceed_idx <- which(misfit > deg2rad(20))
  if (length(exceed_idx) > 0) {
    warning(sprintf(
      "Misfit angle exceeds 20\u00B0 for %d cases. Max misfit: %.2f\u00B0.",
      length(exceed_idx),
      rad2deg(max(misfit))
    ))
  }

  ax <- crossprod(fvec, lvec)
  ang <- (angle(lvec, fvec) - pi / 2) / 2
  list(
    fvec = rotate(fvec, ax, ang),
    lvec = rotate(lvec, ax, -ang),
    misfit = misfit
  )
}

#' @rdname pair_correct
#' @export
correct_pair <- function(x) {
  corr <- misfit_pair(x)
  p <- Plane(corr$fvec)
  l <- Line(corr$lvec)

  if (inherits(x, "fault")) {
    Fault(p[, 1], p[, 2], l[, 1], l[, 2], x[, 5], correction = FALSE)
  } else {
    Pair(p[, 1], p[, 2], l[, 1], l[, 2], correction = FALSE)
  }
}

#' Extract components of fault object
#'
#' [Fault_plane()] extracts the orientation of the fault plane,
#' [Fault_slip()] extracts the orientation of the slip vector, and
#' [Fault_rake()] extracts the rake of the fault, i.e. the angle between fault
#' slip vector and fault strike (measured clockwise from strike!)
#' [Fault_sense()] extracts the fault sense from the rake (1: normal, -1: reverse)
#'
#' @inheritParams Fault_PT
#' @param steps Integer. Either 2, 4, or 8 steps for parsing the fault sense
#'
#' @returns numeric. `"Plane"`, `"Ray"` or angle in degrees, respectively
#' @name Fault_components
#'
#' @examples
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' Fault_plane(f)
#' Fault_slip(f)
#' Fault_rake(f)
#' Fault_sense(f, 2)
#' Fault_sense(f, 4)
#' Fault_sense(f, 8)
NULL

#' @rdname Fault_components
#' @export
Fault_rake <- function(x) {
  stopifnot(is.Fault(x))
  strike <- Line(dd2rhr(x[, 1]), rep(0, nrow(x)))
  slip <- Line(x)
  x[, 5] * angle(strike, slip)
}

#' @rdname Fault_components
#' @export
Pair_plane <- function(x) {
  stopifnot(is.Pair(x))
  # azi <- x[, 3]
  # inc <- x[, 4]
  # sense <- x[, 5]
  # # azi_cor <- ifelse(sense>1, azi + 180, azi)
  # Ray(azi, inc, sense = sense)
  Line(x)
}

#' @rdname Fault_components
#' @export
Fault_slip <- Pair_plane <- function(x) {
  stopifnot(is.Fault(x))
  # azi <- x[, 3]
  # inc <- x[, 4]
  # sense <- x[, 5]
  # # azi_cor <- ifelse(sense>1, azi + 180, azi)
  # Ray(azi, inc, sense = sense)
  Line(x)
}


#' @rdname Fault_components
#' @export
Fault_plane <- Pair_plane <- function(x) {
  stopifnot(is.Pair(x))
  Plane(x)
}

#' @rdname Fault_components
#' @export
Fault_sense <- function(x, steps = 8) {
  rake <- Fault_rake(x) %% 360
  # sign(rake)

  if (steps == 2) {
    # Two groups: N (0–180), R (180–360)
    cut(rake,
      breaks = seq(0, 360, by = 180),
      labels = c("Normal", "Reverse"),
      include.lowest = TRUE,
      right = FALSE
    )
  } else if (steps == 4) {
    # Four groups: NS, ND, RD, RS
    cut(rake,
      breaks = seq(0, 360, by = 90),
      labels = c("Normal-sinistral", "Normal-dextral", "Reverse-dextral", "Reverse-sinistral"),
      include.lowest = TRUE,
      right = FALSE
    )
  } else if (steps == 8) {
    # Eight groups: S, NS, N, ND, D, RD, R, RS
    cut(rake,
      breaks = seq(0, 360, by = 45),
      labels = c("Sinistral", "Normal-sinistral", "Normal", "Normal-dextral", "Dextral", "Reverse-dextral", "Reverse", "Reverse-sinistral"),
      include.lowest = TRUE,
      right = FALSE
    )
  } else {
    stop("steps must be 2, 4, or 8")
  }
}



#' Fault from plane and rake
#'
#' @param p object of class `"Plane"`
#' @param rake Angle (in degrees) between fault strike and lineation.
#' Measured clockwise from the strike, i.e.
#' down is positive; values between 0 and 360&deg; (or −180&deg; and 180&deg;)
#' @param sense Either 1 (for normal fault movement) or -1 (reverse fault movement).
#' Use this only if `rake` values are in range of 0 and 180&deg;.
#' @param ... optional arguments passed to [Fault()]
#'
#' @details Rake is used to describe the direction of fault motion with respect
#' to the strike.
#'
#' Measured clockwise from the strike, down is positive; values between 0 and 360&deg; (or −180&deg; and 180&deg;):
#' - left-lateral strike slip: rake near 0&deg;
#' - right-lateral strike slip: rake near 180&deg;
#' - normal: rake near 90&deg;
#' - reverse/thrust: rake near -90&deg; (270&deg;)
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' fr <- Fault_from_rake(Plane(c(120, 120, 100, 0), c(60, 60, 50, 40)),
#'   rake = c(84.7202, -10, 30, 180)
#' )
#' plot(fr, col = 1:4)
#' legend("topleft", legend = Fault_sense(fr, 8), col = 1:4, pch = 16)
#'
#' fr2 <- Fault_from_rake(Plane(c(90, 90, 90), c(80, 40, 10)),
#'   rake = c(10, 20, 90), sense = c(1, 1, -1)
#' )
#' plot(fr2, col = 1:3)
#' legend("topleft", legend = Fault_sense(fr2, 8), col = 1:3, pch = 16)
Fault_from_rake <- function(p, rake, sense = NULL, ...) {
  stopifnot(is.Plane(p))
  strike <- Ray(dd2rhr(p[, 1]), rep(0, nrow(p)))
  # rake_mod <- rake %% 360
  # rake <- ifelse(rake > 180, rake - 360, rake)
  if (is.null(sense)) {
    qd <- (rake %% 360) <= 180
    sense <- ifelse(qd, 1, -1)
  } else {
    # rake <- sense * (rake %% 180)
  }

  l <- rotate(strike, p, rake)

  Fault(p[, 1], p[, 2], l[, 1], l[, 2], sense = sense, ...)
}

#' Fault from rake and quadrant notation
#'
#' @param p `"Plane"` object
#' @param rake numeric. Rake angle in degrees
#' @param quadrant character. Quadrant of plunge or rake direction
#' @param type character. Either `"plunge"` or `"rake"` for specifying which quadrant convention is used
#' @param sense Either 1 (for normal fault movement) or -1 (reverse fault movement). Only used when `type=="rake"`
#'
#' @details
#' `type=="plunge"` This is the angle measured in the fault plane between the
#' strike given by either right or left-hand rule and the lineation.
#' The angle is recorded in a clockwise sense (looking down upon the fault plane) and has a range
#' from 0 to 180%deg;. The quadrant of plunge indicates the direction
#' of the strike from which the rake angle is measured. The `sense` argument is
#' ignored, as it it implied by the sign of the rake.
#'
#' `type=="rake"` Rake is the **acute** angle measured in the fault plane between the strike of the fault and the
#' lineation . Starting from the strike line, the angle is measured in a sense
#' which is down the dip of the plane. Quadrant of rake indicate the direction of
#' the strike  from which the rake angle is measured, i.e. whether right-hand or left-handrule is followed.
#' Angle ranges from 0 to 90 &deg;. Use `sense` argument to specify the sense of motion.
#'
#' @family parse-orientations
# #' @seealso [azimuth_to_cardinal()] to convert azimuth to cardinal directions, [quadrant2dd()] to define a plane using the strike and dip quadrant notation
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' dip <- c(5, 10, 15, 30, 40, 55, 65, 75, 90)
#' dip_dir <- c(180, 225, 270, 315, 360, 0, 45, 90, 135)
#' rake1 <- c(0, 45, 90, 135, 180, 45, 90, 135, 180)
#' plunge_quadrant <- c("E", "S", "W", "N", "E", "W", "E", "S", "W")
#' Fault_from_rake_quadrant(Plane(dip_dir, dip), rake1, plunge_quadrant, type = "plunge")
#'
#' rake2 <- c(0, 45, 90, 45, 0, 45, 90, 45, 0)
#' rake_quadrant <- c("E", "S", "S", "E", "E", "W", "N", "S", "W")
#' Fault_from_rake_quadrant(Plane(dip_dir, dip), rake2, rake_quadrant, type = "rake")
Fault_from_rake_quadrant <- function(p, rake, quadrant, type = c("plunge", "rake"), sense = NULL) {
  type <- match.arg(type)

  if (type == "plunge") {
    # try right-hand-rule
    f <- Fault_from_rake(p, rake)
    lower_hem <- f[, 4] >= 0
    to_lower <- ifelse(lower_hem, 0, 180)
    azi1 <- f[, 3] + to_lower

    # check if already correct
    res1 <- azimuth_to_cardinal(azi1, n_directions = 16)
    match1 <- sapply(seq_along(rake), function(i) grepl(quadrant[i], res1[i]))
    if (!all(match1)) {
      # do left-handrule if not
      rake2 <- rake + ifelse(match1, 0, -180)
      f <- Fault_from_rake(p, rake2)
    }
  } else {
    dipdir1 <- p[, 1]
    strike1 <- dd2rhr(dipdir1)
    rake_mod <- rake %% 180
    res1 <- azimuth_to_cardinal(strike1, n_directions = 16)
    match1 <- sapply(seq_along(rake_mod), function(i) grepl(quadrant[i], res1[i]))

    rake2 <- 360 + ifelse(match1, 1, -1) * rake
    f <- Fault_from_rake(p, rake2 %% 180, sense = sense)
  }
  return(f)
}


rake2pitch <- function(rake) {
  rake <- abs(rake)
  ifelse(rake > 90, 180 - rake, rake)
}

#' Pitch of a line of a plane
#'
#' @inheritParams Fault_from_rake
#' @param l object of class `"Line"`
#'
#' @returns Pitch angle in degrees
#' @export
#'
#' @examples
#' Plane_pitch(Plane(90, 70), Line(175, 13)) # 13.9
Plane_pitch <- function(p, l) {
  Fault_rake(Fault(p[, 1], p[, 2], l[, 1], l[, 2], 1)) |> rake2pitch()
}

#' Apparent dip direction
#'
#' @inheritParams Fault_from_rake
#' @param apparent_dip angle in degrees
#'
#' @returns Azimuth in degrees
#' @export
#'
#' @examples
#' apparent_dip_direction(Plane(63, 45), 10)
apparent_dip_direction <- function(p, apparent_dip) {
  res <- Fault_from_rake(p, rake = -apparent_dip)
  res[, 3]
}

#' Apparent dip direction
#'
#' @param a1,a2 `"Line"` objects containing the apparent dips and dip directions
#'
#' @returns `"Plane"` object
#' @export
#'
#' @examples
#' a1 <- Line(45, 22)
#' a2 <- Line(352, 10)
#' res <- Plane_from_apparent_dips(a1, a2)
#'
#' stereoplot()
#' points(rbind(a1, a2))
#' lines(res, lty = 2)
Plane_from_apparent_dips <- function(a1, a2) {
  res <- crossprod(a1, a2)
  Plane(res)
}



#' Dip-slip Kinematics from Strike-Slip Faults
#'
#' Dip-slip sense (1 for normal, -1 for reverse) when only strike-slip kinematics are known
#'
#' @param x `"Pair"` object(s) of the plane(s) and line(s) with unknown dip-slip offset
#' @param left logical. `TRUE` if `x` is a left-lateral (sinistral) fault, and
#' `FALSE` if `x` is a right-lateral (dextral) fault.
#' Must have same length as rows in `x`
#'
#' @returns numeric. 1 if normal, -1 if reverse offset
#' @family parse-orientations
#' @seealso [strikeslip_kinematics()]
#'
#' @export
#'
#' @examples
#' # Sinistral fault
#' sense_from_strikeslip(Pair(120, 89, 30, 5), left = TRUE) # 1: normal offset
#'
#' # Dextral fault
#' sense_from_strikeslip(Pair(120, 89, 30, 5), left = FALSE) # -1: reverse offset
sense_from_strikeslip <- function(x, left) {
  stopifnot(is.Pair(x), is.logical(left), length(left) == nrow(x))

  strike <- Line(dd2rhr(x[, 1]), rep(0, nrow(x)))
  slip <- Line(x)
  rake <- angle(strike, slip)

  # if sinistral fault:
  sense <- ifelse(rake < 90, 1, -1)

  # opposite sign if dextral:
  ifelse(left, sense, -sense)
}


#' Strike-slip Kinematics
#'
#' Returns the strike-slip kinematics of a fault
#'
#' @param x `"Fault"` object(s)
#'
#' @returns character. `"left"` - left-lateral (sinistral) offset, `"right"` - right-lateral (dextral) offset
#' @export
#' @family parse-orientations
#' @seealso [sense_from_strikeslip()]
#'
#' @examples
#' strikeslip_kinematics(Fault(120, 50, 104, 49, sense = -1))
strikeslip_kinematics <- function(x) {
  stopifnot(is.Fault(x))
  normal <- x[, 5] == 1

  strike <- Line(dd2rhr(x[, 1]), rep(0, nrow(x)))
  slip <- Line(x)
  rake <- angle(strike, slip)

  # if normal fault
  sinistral <- rake < 90

  # dextral if reverse fault
  sinistral <- ifelse(normal, sinistral, !sinistral)

  # return string
  ifelse(sinistral, "left", "right")
}
