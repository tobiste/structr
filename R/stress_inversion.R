#' Stress Inversion for Fault-Slip Data
#' 
#' Determines the orientation of the principal stresses from fault slip data using the Michael (1984) method.
#' Confidence intervals are estimated by bootstrapping. 
#' This inversion is simplified by the assumption that the magnitude of the 
#' tangential traction on the various fault planes, at the time of rupture, is similar.
#' 
#' @param x `"Fault"` object
#' @param boot integer. Number of bootstrap samples (10 by default)
#' @param conf.level numeric. Confidence level of the interval (0.95 by default)
#' 
#' @returns list
#'  \describe{
#'  \item{`stress_tensor`}{matrix. Best-fit devitoric stress tensor}
#'  \item{`principal_axes`}{`"Line"` obects. Orientation of the principal stress axes}
#'  \item{`principal_axes_conf`}{list containg the confidence ellipses for the 3 principal stress vectors. See [confidence_ellipse()] for details.}
#'  \item{`principal_vals`}{numeric. The proportional magnitudes of the principal stress axes given by the eigenvalues of the stress tensor: \eqn{\sigma_1}, \eqn{\sigma_2}, and \eqn{\sigma_3}}
#'  \item{`principal_vals_conf`}{3-colum vector containing the lower and upper margins of the confidence interval of the principal vals}
#'  \item{`R`}{numeric. Stress shape ratio after Gephart & Forsyth (1984): \eqn{R = (\sigma_1 - \sigma_2)/(\sigma_1 - \sigma_3)}. Values ranging from 0 to 1, with 0 being
#' \eqn{\sigma_1 = \sigma_2} and 1 being \eqn{\sigma_2 = \sigma_3}.}
#'  \item{`R_conf`}{Confidence interval for `R`}
#'  \item{`phi`}{numeric. Stress shape ratio after Angelier (1979): \eqn{\Phi = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)}. Values range between 0 (\eqn{\sigma_2 = \sigma_3}) and 1 (\eqn{\sigma_2 = \sigma_1}).}
#'  \item{`phi_conf`}{Confidence interval for `phi`}
#'  \item{`bott`}{numeric. Stress shape ratio after Bott (1959): \eqn{\R = (\sigma_3 - \sigma_1)/(\sigma_2 - \sigma_1)}. Values range between \eqn{-\infty} and \eqn{+\infty}.}
#'  \item{`bott_conf`}{Confidence interval for `bott`}
#'  \item{`beta`}{numeric. Average angle between the tangential traction predicted by the best stress tensor  and the slip vector on each plane. Should be close to 0.}  
#'  \item{`sigma_s`}{numeric. Average resolved shear stress on each plane. Should be close to 1.}
#'  \item{`fault_data`}{`data.frame` containing the beta angles, the angles between sigma 1 and the plane normal, 
#'  the resolved shear and normal stresses, and the slip and dilation tendency on each plane.}
#'  }
#' 
#' @references Michael, A. J. (1984). Determination of stress from slip data: Faults and folds. Journal of Geophysical Research: Solid Earth, 89(B13), 11517–11526. https://doi.org/10.1029/JB089iB13p11517
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
#' 
#' 
#' @export
#' 
#' @importFrom stats t.test
#' 
#' @seealso [SH()] to calculate the azimuth of the maximum horizontal stress; [Fault_PT()] for a simple P-T stress analysis.
#' 
#' @examples
#' # Use Angelier examples:
#' res_TYM <- slip_inversion(angelier1990$TYM, boot = 10)
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
#' res_AVB <- slip_inversion(angelier1990$AVB)
#' stereoplot(title = "Agia Varvara, Crete, Greece", guides = FALSE)
#' fault_plot(angelier1990$AVB, col = "gray80")
#' stereo_confidence(res_AVB$principal_axes_conf$sigma1, col = 2)
#' stereo_confidence(res_AVB$principal_axes_conf$sigma2, col = 3)
#' stereo_confidence(res_AVB$principal_axes_conf$sigma3, col = 4)
#' text(res_AVB$principal_axes, label = rownames(res_AVB$principal_axes), col = 2:4, adj = -.25)
#' legend("topleft", col = 2:4, legend = rownames(res_AVB$principal_axes), pch = 16)
slip_inversion <-  function(x, boot = 10L, conf.level = 0.95){
  best.fit <- slip_inversion0(x)
  fault_df <- best.fit$fault_data
  nx <- nrow(x)
  
  # bootstrap results
  boot_results <- lapply(1:boot, function(i) {
    idx <- sample.int(nx, replace = TRUE)
    x_sample <- x[idx, ]
    slip_inversion0(x_sample)
  })
  
  # calculate confidence intervals from bootstrap results
  sigma_vec1 <- do.call(rbind, lapply(boot_results, function(x){
    x$principal_axes[1, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, res = 100)
  
  sigma_vec2 <- do.call(rbind, lapply(boot_results, function(x){
    x$principal_axes[2, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, res = 100)
  
  sigma_vec3 <- do.call(rbind, lapply(boot_results, function(x){
    x$principal_axes[3, ]
  })) |> confidence_ellipse(alpha = 1 - conf.level, res = 100)
  
  R_boot <- vapply(boot_results, function(x){x$R}, FUN.VALUE = numeric(1)) |> 
    stats::t.test(conf.level = conf.level)
  
  phi_boot <- vapply(boot_results, function(x){x$phi}, FUN.VALUE = numeric(1)) |> 
    stats::t.test(conf.level = conf.level)
  
  bott_boot <- vapply(boot_results, function(x){x$bott}, FUN.VALUE = numeric(1)) |> 
    stats::t.test(conf.level = conf.level)
  
  sigma_boot0 <- vapply(boot_results, function(x){x$principal_vals}, FUN.VALUE = numeric(3)) |> t()
  sigma_boot <- sapply(1:3, function(col){
    sigma_boot_col <- stats::t.test(sigma_boot0[, col], conf.level = conf.level)
    sigma_boot_col$conf.int
  })
  colnames(sigma_boot) <- names(best.fit$principal_vals)
  
  list(
    stress_tensor = best.fit$stress_tensor,
    principal_axes = best.fit$principal_axes,
    principal_axes_conf = list(sigma1 = sigma_vec1, sigma2 = sigma_vec2, sigma3 = sigma_vec3),
    principal_vals = best.fit$principal_vals,
    principal_vals_conf = sigma_boot,
    R = best.fit$R,
    R_conf = R_boot$conf.int,
    phi = best.fit$phi,
    phi_conf = phi_boot$conf.int,
    bott = best.fit$bott,
    bott_conf = bott_boot$conf.int,
    beta = best.fit$beta,
    sigma_s = best.fit$sigma_s,
    fault_data = fault_df
  )
}

slip_inversion0 <- function(x){
  tau <- linear_stress_inversion(x)
  # tau0 <- tau / sqrt(sum(tau^2)) # normalize Frobenius norm
  
  # Eigen decomposition of stress tensor
  eig <- eigen(tau)
  # sigma_vals <- sort(eig$values, decreasing  = TRUE)
  sigma_vals <- eig$values

  # stress ratios:
  R <- (sigma_vals[1] - sigma_vals[2]) / (sigma_vals[1] - sigma_vals[3]) # Gephart & Forsyth 1984
  phi <- (sigma_vals[2] - sigma_vals[3]) / (sigma_vals[1] - sigma_vals[3]) # Angelier 1979
  shape_ratio_bott <- (sigma_vals[3] - sigma_vals[1]) / (sigma_vals[2] - sigma_vals[1]) # Bott, Simon-Gomez

  principal_axes <- t(eig$vectors) |> as.Vec3() |> Line() # sigma1, sigma2, sigma3
  names(sigma_vals) <- rownames(principal_axes) <- c("sigma1", 'sigma2', "sigma3")
  
  
  # Angles between the tangential traction predicted by the best stress tensor and the slip vector on each plane
  betas <- sapply(1:nrow(x), function(i){
    int <- crossprod(Plane(principal_axes[2, ]), Plane(x[i, ])) |> Line()
    angle(int, Line(x[i, ]))
  }) #|> as.vector()
  betas <- ifelse(betas > 90, 180 - betas, betas)
  beta_mean <- tectonicr::circular_mean(betas)
  
  
  # Resolved shear stress on plane
  theta <- sapply(1:nrow(x), function(i){
    angle(Plane(x[i, ]), principal_axes[1, ])
  })
  
  sigma_s <- shear_stress(sigma_vals[1], sigma_vals[3], theta)
  sigma_n <- normal_stress(sigma_vals[1], sigma_vals[3], theta)
  slip_tend <- slip_tendency(sigma_s, sigma_n)
  dilat_tend <- dilatation_tendency(sigma_vals[1], sigma_vals[3], sigma_n)
  
  sigma_s_mean <- mean(sigma_s)

  list(
    stress_tensor = tau,
    principal_axes = principal_axes,
    principal_vals = sigma_vals,
    R = R,
    phi = phi,
    bott = shape_ratio_bott,
    beta = beta_mean,
    sigma_s = sigma_s_mean,
    fault_data = data.frame(beta=betas, theta=theta, sigma_s=sigma_s, sigma_n = sigma_n, slip_tendency = slip_tend, dilational_tendency = dilat_tend)
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
  #least_square_solv <- lm(as.vector(t(s)) ~ 0 + A)
  #stress_vector <- least_square_solv$coefficients
  
  # trace of stress tensor tau is usually assumed to be zero (Michael 1984):
  # trace(tau)  = sigma1 + sigma2 + sigma3 = 0 (Eq. 2)
  
  # Taken in combination with the constraint that the isotropic stress is zero, i.e. s33 = -(s11 + s22))
  stress_vector[6] <- -(stress_vector[1] + stress_vector[4])
  names(stress_vector) <- c('11', '12', '13', '22', '23', '33')
  
  stress_tensor <- matrix(c(
    stress_vector['11'], stress_vector['12'], stress_vector['13'],
    stress_vector['12'], stress_vector['22'], stress_vector['23'],
    stress_vector['13'], stress_vector['23'], stress_vector['33']
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
