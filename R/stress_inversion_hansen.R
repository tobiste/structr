#' 9D Direct Inversion for Fault Slip Including Vorticity
#'
#' Direct inversion of stress, strain or strain rate including vorticity using
#' 9D parameter space using the method by Hansen (2013). It can be applied
#' regardless whether the dynamic or the kinematic hypothesis is adopted;
#' it can handle datasets representing two to seven degrees of freedom; and it
#' is not dependent on the correct assessment of slip sense.  
#' If no vorticity is involved, the inversion can be done by using a 6-dimensional 
#' parameter space only (`type = '6d'`).
#'
#' @param x object of class `"Pair"` or `"Fault"` with at least 7 rows.
#' @param flip logical. Flip if you want to have the negative stress tensor, i.e.
#' sigma 1 and 3 will be flipped.
#' @param type character. Inversion method, either `"9d"` (the default) for using the 
#' 9-dimensional or `"6d"` for the 6-dimensional parameter space.
#' @inheritParams slip_inversion_michael
#'
#' @returns list. See [slip_inversion_michael()] for output description. 
#' If `type == '9d`, additional outputs are
#' the vorticity axis (`"vorticity_axis"`, a `Vec3` object) and the magnitude of
#' vorticity (`"vorticity_mag"`, a numeric).
#' 
#' @details
#' ## Pole to the M-plane
#' \deqn{\mathbf{b} = \mathbf{n} \times \mathbf{v}}
#' where \eqn{\mathbf{n}} is the upward unit normal to the fault plane and \mathbf{v}
#'  is the unit slip vector.
#' 
#' ## 9D f-poles
#' \deqn{\mathbf{f}_{nr} = \left[ b_1 n_1, b_1 n_2, b_1 n_3, b_2 n_1, b_2 n_2, b_2 n_3, b_3 n_1, b_3 n_2, b_3 n_3\right]}
#' \deqn{\hat{\mathbf{f}} = \frac{\mathbf{f}_{nr}}{|\mathbf{f}_{nr}|}}
#' 
#' ## Second moment tensor
#' \deqn{\hat{M} = \sum_{r = 1}^{N} \hat{\mathbf{f}}  \otimes \hat{\mathbf{f}} }
#' 
#' ## Inverted slip tensor
#' The 9D stress vector \eqn{\hat{s}} is the eigenvector of \eqn{\hat{M}} 
#' corresponding to the second-lowest eigenvalue, reshaped into the asymmetric 
#' inverted slip tensor:
#' 
#' \deqn{\hat{\dot{T}} = \begin{pmatrix}
#'    \hat{s}_1 & \hat{s}_2 & \hat{s}_3 \\
#'    \hat{s}_4 & \hat{s}_5 & \hat{s}_6 \\
#'    \hat{s}_7 & \hat{s}_8 & \hat{s}_9 \\
#' \end{pmatrix}
#'  }
#' 
#' ## Symmetric and antisymmetric decomposition
#' \deqn{\hat{\dot{T}}_S = \frac{\hat{\dot{T}} + \hat{\dot{T}}^{\top}}{2} }
#' 
#' \deqn{\hat{\dot{T}}_A = \frac{\hat{\dot{T}} - \hat{\dot{T}}^{\top}}{2} }
#' 
#' ## Principal axes and shape ratio
#' 
#' Eigen-decompose \eqn{\hat{\dot{T}}_S}, sort eigenvalues descending \eqn{\lambda_1 \geq \lambda_2 \geq \lambda_3}. 
#' The eigenvectors give the principal stress axes \eqn{\mathbf{s}_1}, \eqn{\mathbf{s}_2}, \eqn{\mathbf{s}_3}. 
#' The shape ratio is:
#' 
#' \deqn{\phi = \frac{\lambda_2 - \lambda_3}{\lambda_1 - \lambda_3}}
#' 
#' ## Reduced symmetric tensor 
#' 
#' \deqn{\mathbf{T}_2 = \mathbf{V} 
#' \begin{pmatrix}
#'    1 & 0 & 0 \\
#'    0 & \phi & 0 \\
#'    0 & 0 & 0
#' \end{pmatrix}
#' \mathbf{V}^{\top}
#' }
#'  where \eqn{\mathbf{V} = [\mathbf{s}_1\ \mathbf{s}_2\ \mathbf{s}_3]} has the 
#'  eigenvectors as columns.
#'  
#'  ## Normalise the antisymmetric part
#'  
#'  \deqn{\hat{T}_A = \hat{\dot{T}}_A \odot \frac{\mathbf{T}_S}{\hat{\dot{T}}_S} }
#'  
#'  where \eqn{\odot} denotes element-wise multiplication and division.
#'  
#'  ## Vorticity axis and magnitude 
#'  
#'  The axial vector \eqn{\hat{T}_A} is
#'  \deqn{\overrightarrow{\omega} = \begin{pmatrix}
#'    \hat{T}_{A,32} \\
#'    \hat{T}_{A,13} \\
#'    \hat{T}_{A,21}
#'  \end{pmatrix}
#'  }
#'  
#'  The unit vorticity axis in geographic coordinates:
#'  \deqn{\mathbf{u}_{xyz} = \frac{\overrightarrow{\omega}}{| \overrightarrow{\omega} |}}
#'  
#'  The vorticity magnitude:
#'  \deqn{|\omega| = 2 | \overrightarrow{\omega} |}
#'
#' @references
#' Hansen, J. A. (2013). Direct inversion of stress, strain or strain rate
#' including vorticity: A linear method of homogenous fault-slip data inversion
#' independent of adopted hypothesis. Journal of Structural Geology, 51, 3–13.
#' https://doi.org/10.1016/j.jsg.2013.03.014
#'
#' @returns list
#' @export
#'
#' @family stress-inversion
#'
#' @examples
#' # Osmundsen et al. 2010 dataset
#' ## 9D solution
#' res <- slip_inversion_hansen(osmundsen2010, flip = TRUE)
#'
#' phi_val <- round(res$stress_shape$phi, 2)
#' rup_val <- round(res$misfit$rup, 2)
#' w_val <- round(res$vorticity_mag, 2)
#'
#' stereoplot(title = "9D: Lofoten / Northern Norway\n(Osmundsen et al. 2010)", guides = FALSE)
#' stereo_shmax(res$SHmax)
#' points(Plane(osmundsen2010), col = assign_col(res$misfit$rup), pch = 0, cex = 0.5)
#' points(Line(osmundsen2010), col = assign_col(res$misfit$rup), pch = 16, cex = 0.5)
#' points(res$principal_axes, col = 2:4, pch = 16, cex = 2)
#' text(res$principal_axes, labels = rownames(res$principal_axes), col = 2:4, adj = -.5)
#' points(res$vorticity_axis, col = 5, pch = 17, cex = 2)
#' text(res$vorticity_axis, labels = bquote(omega), col = 5, adj = -.5)
#'
#' title(sub = bquote(Phi == .(phi_val) ~ "|" ~ bar("RUP") == .(rup_val) * "%" ~
#'   "|" ~ omega == .(w_val)))
#'
#' ## 6D inversion
#' res6 <- slip_inversion_hansen(osmundsen2010, flip = TRUE, type = "6d")
#'
#' phi6_val <- round(res6$stress_shape$phi, 2)
#' rup6_val <- round(res6$misfit$rup, 2)
#'
#' stereoplot(title = "6D: Lofoten / Northern Norway\n(Osmundsen et al. 2010)", guides = FALSE)
#' stereo_shmax(res6$SHmax)
#' points(Plane(osmundsen2010), col = assign_col(res6$misfit$rup), pch = 0, cex = 0.5)
#' points(Line(osmundsen2010), col = assign_col(res6$misfit$rup), pch = 16, cex = 0.5)
#' points(res6$principal_axes, col = 2:4, pch = 16, cex = 2)
#' text(res6$principal_axes, labels = rownames(res6$principal_axes), col = 2:4, adj = -.5)
#'
#' title(sub = bquote(Phi == .(phi6_val) ~ "|" ~ bar("RUP") == .(rup6_val) * "%"))
#'
#' # Angelier 1990 dataset
#' nx <- length(angelier1990)
#' par(mfrow = c(1, nx))
#'
#' invisible(lapply(seq_len(nx), function(i) {
#'   # inversion
#'   x <- angelier1990[[i]]
#'   res <- slip_inversion_hansen(x, type = "6d")
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
slip_inversion_hansen <- function(x, flip = FALSE, friction = 0.6, type = c("9d", "6d")) {
  tsign <- if (flip) -1 else 1
  type <- match.arg(type)

  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x))
  slips <- if (is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips)

  n <- nrow(normals)
  if (n < 7L) {
    warning(
      "Fewer than 7 fault-slip measurements: solution may be ",
      "underdetermined."
    )
  }

  # 1. Compute poles to M-planes and f-poles  (Eqs. 14, 19)
  # n_vec = fault plane normal  (n in paper)
  # b_vec = pole to M-plane     (b in paper) = n x v  (Fig. 1)
  # f-pole (bf, Eq. 19): bf = [b1n1, b1n2, b1n3, b2n1, ..., b3n3] / |.|
  b <- crossprod(normals, slips)

  n1 <- normals[, 1]
  n2 <- normals[, 2]
  n3 <- normals[, 3]

  b1 <- b[, 1]
  b2 <- b[, 2]
  b3 <- b[, 3]

  if (type == "9d") {
    f <- cbind(
      b1 * n1, b1 * n2, b1 * n3,
      b2 * n1, b2 * n2, b2 * n3,
      b3 * n1, b3 * n2, b3 * n3
    )
    f <- f / sqrt(rowSums(f * f))

    # 2. Second moment tensor  (Eq. 21)
    M <- crossprod(f)

    # 3. Eigen-decomposition; second eigenvector = best estimate of the slip vector
    #    (Eq. 20)
    eM <- eigen(M, symmetric = TRUE)
    MSort <- order(Re(eM$values)) # ascending
    val_M <- Re(eM$values)[MSort]
    vec_M <- t(Re(eM$vectors))[MSort, ] # rows = eigenvectors

    s <- vec_M[2, ] # second eigenvector in ascending rank

    # 4. Build inverted slip tensor; symmetric + antisymmetric decomposition
    #    (Eqs. 8–10; algorithm steps 6–7)
    T_mat <- matrix(s, 3, 3, byrow = TRUE) * tsign
    Ts <- (T_mat + t(T_mat)) / 2 # b_T_hat_S
    Ta <- (T_mat - t(T_mat)) / 2 # b_T_hat_A


    # 5. Principal axes and shape ratio from symmetric part  (Eqs. 2–5)
    eTs <- eigen(Ts, symmetric = TRUE)
    val_Ts <- eTs$values
    vec_Ts <- eTs$vectors

    phi <- (val_Ts[2] - val_Ts[3]) / (val_Ts[1] - val_Ts[3])

    # Reconstruct reduced symmetric tensor T_S = V D V^T  (Eq. 5)
    # V has eigenvectors as columns = t(vec_Ts)
    # Norm <- t(vec_Ts) %*% diag(c(nval_Ts[1], nval_Ts[2], 0)) %*% vec_Ts

    # Because the diagonal is (1, phi, 0), the Eq. is actually equivalent to:
    v1 <- vec_Ts[, 1]
    v2 <- vec_Ts[, 2]
    Norm <- tcrossprod(v1) + phi * tcrossprod(v2)

    # 6. Vorticity axis and magnitude  (Eqs. 12–13, 22)
    # Normalise antisymmetric part: b_T_A = b_T_hat_A * (T_S / b_T_hat_S)
    ratio <- Norm / Ts
    ratio[Ts == 0 & Norm == 0] <- 0 # guard 0/0 -> 0
    Ta <- Ta * ratio

    # Axial vector of T_A (Eq. 10)
    w <- as.Vec3(c(Ta[3, 2], Ta[1, 3], Ta[2, 1]))
    w_xyz <- vnorm(w)
    w_mag <- -2 * vlength(w)

    # Enforce downward-pointing convention
    # if (nw[3] >= 0) nw <- -nw else wlen <- -wlen
  } else {
    # --- 6D f-pole (Eq. 15): f = [b1n1, b2n2, b3n3, b1n2+b2n1, ...] ---
    f6 <- cbind(
      b1 * n1,
      b2 * n2,
      b3 * n3,
      b1 * n2 + b2 * n1,
      b2 * n3 + b3 * n2,
      b1 * n3 + b3 * n1
    )
    f6 <- f6 / sqrt(rowSums(f6 * f6))

    M6 <- crossprod(f6)

    eM6 <- eigen(M6, symmetric = TRUE)
    M6Sort <- order(eM6$values)
    val_M6 <- eM6$values[M6Sort]
    vec_M6 <- t(eM6$vectors)[M6Sort, ]

    s6 <- vec_M6[2, ] # second eigenvector in ascending rank (Eq. 16)


    # 6D: symmetric by construction (Eq. 16)
    T6_mat <- matrix(
      c(
        s6[1], s6[4], s6[6],
        s6[4], s6[2], s6[5],
        s6[6], s6[5], s6[3]
      ),
      nrow = 3,
      byrow = TRUE
    )
    Ts <- T6_mat * tsign

    eT6 <- eigen(Ts, symmetric = TRUE)
    val_Ts <- eT6$values
    vec_Ts <- eT6$vectors
  }

  principal_axes <- as.Vec3(t(vec_Ts))
  rownames(principal_axes) <- names(val_Ts) <- c("sigma1", "sigma2", "sigma3")

  # 8.
  stress_shape <- stress_shape(Ts)

  # --- Step 4: Per-fault diagnostics ---
  misfit <- slip_inversion_misfit(Ts, x)

  # Angle between slip planes and sigma 1
  theta <- vapply(seq_len(n), function(i) {
    angle(Plane(x[i, ]), principal_axes[1, ])
  }, numeric(1))

  # Theoretically resolved shear stress on plane
  sigma_s_mean <- mean(abs(shear_stress(val_Ts[1], val_Ts[3], theta)))

  shearnorm <- tau2shearnorm(Ts, x, friction = friction)
  tendency <- tau2tendency(Ts, x, friction = friction)

  SHmax <- tryCatch(
    expr = SH_from_tensor(eigen(Ts, symmetric = TRUE)$vectors),
    error = function(e) {
      SH(principal_axes[1, ], principal_axes[2, ], principal_axes[3, ], R = stress_shape$R)
    }
  )

  pfaults <- principal_fault(principal_axes[1, ], principal_axes[3, ], friction)


  # 9 results

  res <- list(
    stress_tensor = as.ellipsoid(Ts),
    principal_axes = Line(principal_axes),
    principal_vals = val_Ts,
    principal_faults = pfaults,
    stress_shape = stress_shape,
    tau_mean = sigma_s_mean,
    stress_components = cbind(shearnorm, tendency),
    misfit = misfit,
    SHmax = SHmax
  )

  if (type == "9d") {
    append(res, list(
      vorticity_mag = w_mag,
      vorticity_axis = w,
      method = "hansen (9d)"
    ))
  } else {
    append(res, list(method = "hansen (6d)"))
  }
}

#' 9D Direct Inversion for Fault Slip Including Vorticity with Confidence Intervals
#'
#' Fault-slip inversion method after Hansen (2013) with bootstrapped confidence intervals
#'
#' @inheritParams slip_inversion_hansen
#' @inheritParams slip_inversion_michael
#'
#' @returns See [slip_inversion_hansen()] and [slip_inversion_michael()]
#' @export
#' @family stress-inversion
#'
#' @examples
#' set.seed(20250411)
#' res <- slip_inversion_hansen_boot(osmundsen2010, n_iter = 100, n = 1000, res = 100)
#'
#' # some stress shape
#' phi_val <- round(res$phi_CI, 2)
#'
#' # vorticity
#' w_val <- round(res$vorticity_mag_CI, 2)
#'
#' # Plot the faults
#' plot(osmundsen2010, col = "grey", lwd = 0.1, cex = 0.5)
#' stereo_confidence(res$principal_axes_CI$sigma1, col = 2)
#' stereo_confidence(res$principal_axes_CI$sigma2, col = 3)
#' stereo_confidence(res$principal_axes_CI$sigma3, col = 4)
#' stereo_confidence(res$vorticity_axis_CI, col = 5)
#' text(res$principal_axes, label = rownames(res$principal_axes), col = 2:4, adj = -.25)
#' text(res$vorticity_axis, labels = bquote(omega), col = 5, adj = -.5)
#' title(
#'   main = "Lofoten / Northern Norway\n(Osmundsen et al. 2010)",
#'   sub = bquote(atop(
#'     varphi ~ "(95% CI)" == "[" * .(phi_val[1]) * "," ~ .(phi_val[2]) * "]",
#'     ~omega ~ "(95%)" == "[" * .(w_val[1]) * "," ~ .(w_val[2]) * "]"
#'   ))
#' )
slip_inversion_hansen_boot <- function(
  x,
  friction = 0.6,
  flip = FALSE,
  type = c("9d", "6d"),
  n_iter = 100,
  conf.level = 0.95,
  ...
) {
  type <- match.arg(type)
  best.fit <- slip_inversion_hansen(x, friction = friction, flip = flip, type = type)

  if (n_iter == 0) {
    return(best.fit)
  } else {
    normals <- unclass(Vec3(Plane(x)))
    slips <- if (is.Fault(x)) Ray(x) else Line(x)
    slips <- unclass(Vec3(slips))

    nx <- nrow(normals)

    tsign <- if (flip) -1 else 1

    # bootstrap results
    res_boot <- future.apply::future_lapply(seq_len(n_iter), function(i) {
      idx <- sample.int(nx, replace = TRUE)
      x_sample <- x[idx, ]
      res <- slip_inversion_hansen(x[idx, ], friction = friction, flip = flip, type = type)
      tau <- res$stress_tensor
      w <- res$vorticity_axis
      list(tau = tau, w = w)
    }, future.seed = TRUE)

    tau_boot <- lapply(res_boot, function(x) x$tau)

    princ_boot <- lapply(tau_boot, tau2stress)


    # calculate confidence intervals from bootstrap results
    sigma_vec1 <- do.call(rbind, lapply(princ_boot, function(x) {
      x$principal_axes[1, ]
    })) |>
      confidence_ellipse(alpha = 1 - conf.level, ...)

    sigma_vec2 <- do.call(rbind, lapply(princ_boot, function(x) {
      x$principal_axes[2, ]
    })) |>
      confidence_ellipse(alpha = 1 - conf.level, ...)

    sigma_vec3 <- do.call(rbind, lapply(princ_boot, function(x) {
      x$principal_axes[3, ]
    })) |>
      confidence_ellipse(alpha = 1 - conf.level, ...)


    # Vorticity
    if (type == "9d") {
      vort_boot <- lapply(res_boot, function(x) x$w)
      vorticity_mag <- do.call(rbind, lapply(vort_boot, function(x) {
        -2 * vlength(x)
      })) |>
        stats::t.test(conf.level = conf.level)

      vort_vec3 <- do.call(rbind, vort_boot) |>
        confidence_ellipse(alpha = 1 - conf.level, ...)
    }

    params_boot <- lapply(tau_boot, stress_shape)

    R_boot <- vapply(params_boot, function(x) {
      x$R
    }, FUN.VALUE = numeric(1)) |>
      stats::t.test(conf.level = conf.level)

    phi_boot <- 1 - rev(R_boot$conf.int)


    bott_boot <- vapply(params_boot, function(x) {
      x$bott
    }, FUN.VALUE = numeric(1)) |>
      stats::t.test(conf.level = conf.level)

    sigma_boot0 <- vapply(princ_boot, function(x) {
      x$sigma_vals
    }, FUN.VALUE = numeric(3)) |>
      t()
    sigma_boot <- sapply(1:3, function(col) {
      sigma_boot_col <- stats::t.test(sigma_boot0[, col], conf.level = conf.level)
      sigma_boot_col$conf.int
    })
    colnames(sigma_boot) <- names(best.fit$principal_vals)
    attr(sigma_boot, "conf.level") <- conf.level


    alpha_CI <- stats::t.test(best.fit$misfit$alpha, conf.level = conf.level)
    rup_CI <- stats::t.test(best.fit$misfit$rup, conf.level = conf.level)

    SHmax_CI <- vapply(tau_boot, function(x) {
      tryCatch(
        expr = SH_from_tensor(x),
        error = function(e) {
          phi <- stress_shape(x)$phi
          principal_axes <- tau2stress(x)$principal_axes
          SH(principal_axes[1, ], principal_axes[2, ], principal_axes[3, ], R = phi)
        }
      )
    }, FUN.VALUE = numeric(1)) |>
      stats::t.test(conf.level = conf.level)

    res <- append(best.fit, list(
      principal_axes_CI = list(sigma1 = sigma_vec1, sigma2 = sigma_vec2, sigma3 = sigma_vec3),
      principal_vals_CI = sigma_boot,
      SHmax_CI = SHmax_CI$conf.int,
      R_CI = R_boot$conf.int,
      phi_CI = phi_boot,
      bott_CI = bott_boot$conf.int,
      alpha_CI = alpha_CI$conf.int,
      rup_CI = rup_CI$conf.int
    ))

    if (type == "9d") {
      append(res, list(
        vorticity_mag_CI = vorticity_mag$conf.int,
        vorticity_axis_CI = vort_vec3
      ))
    } else {
      res
    }
  }
}
