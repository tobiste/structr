#' Fault slip stress inversion
#' 
#' Michael (1984) method
#' 
#' @param x `"Fault"` object
#' @param friction numeric. Value(s) for coefficient of friction
#' 
#' @returns list
#' 
#' @references Michael, A. J. (1984). Determination of stress from slip data: Faults and folds. Journal of Geophysical Research: Solid Earth, 89(B13), 11517–11526. https://doi.org/10.1029/JB089iB13p11517
#' @export
#' 
#' @examples
#' stress_inversion(angelier1990$TYM)
#' stress_inversion(angelier1990$AVB)
stress_inversion <- function(x, friction = 0.6){
  tau <- linear_stress_inversion(x)
  # tau0 <- tau / sqrt(sum(tau^2)) # normalize Frobenius norm
  
  # Eigen decomposition of stress tensor
  eig <- eigen(tau)
  #sigma_vals <- sort(eig$values, decreasing  = TRUE)
  sigma_vals <- eig$values

  # stress ratios:
  R <- (sigma_vals[1] - sigma_vals[2]) / (sigma_vals[1] - sigma_vals[3]) # Gephart & Forsyth 1984
  phi <- (sigma_vals[2] - sigma_vals[3]) / (sigma_vals[1] - sigma_vals[3]) # Angelier 1979
  shape_ratio_bott <- (sigma_vals[3] - sigma_vals[1]) / (sigma_vals[2] - sigma_vals[1]) # Bott, Simon-Gomez

  #maybe transpose?
  principal_axes <- t(eig$vectors) |> as.Vec3() |> Line() # sigma1, sigma2, sigma3
  names(sigma_vals) <- rownames(principal_axes) <- c("sigma1", 'sigma2', "sigma3")
  
  # Principal stress directions
  # sigma_vec1 <- principal_axes[1, ]
  # sigma_vec2 <- principal_axes[2, ]
  # sigma_vec3 <- principal_axes[3, ]
  
  # p <- Plane(x)
  # 
  # mean_instability <- numeric(length(friction))
  # 
  # instabilites <- lapply(friction, function(i){
  #   fault_instability_criterion(x, R = shape_ratio_gephart, mu = i)
  # })
  # mean_instability <- sapply(instabilites, mean)
  
  
  list(
    stress_tensor = tau0,
    principal_axes = principal_axes,
    sigma_vals = sigma_vals,
    R = R,
    phi = phi,
    bott = shape_ratio_bott
    #friction = friction,
    #instability = instabilites,
    #mean_instability = mean_instability
  )
}

linear_stress_inversion <- function(fault) {
  m <- nrow(fault)
  
  n <- Plane(fault) |>
    Vec3() |>
    unclass() # plane normal
  s <- Line(fault) |>
    Vec3() |>
    unclass() # slip vector
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
  
  # not sure if this is the correct configuration:
  # stress_tensor <- matrix(c(
  #   stress_vector['11'], stress_vector['12'], stress_vector['13'], # sigma 1
  #   0, stress_vector['22'], stress_vector['23'], # sigma 2
  #   0, 0, stress_vector['33'] # sigma 3
  # ), nrow = 3, byrow = TRUE)
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
