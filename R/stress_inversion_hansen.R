# FSinvert - Fault-slip data inversion
# Translated from Python by John Are Hansen (GPL-3) to R.
#
# Input:  tab/space-delimited text file with columns [Strike, Dip, Trend, Plunge]
#         (no header). Strike in RHR convention, all values in degrees.
#         
#         Note that the dataset represents the slip vector, not necessarily correct slip sense. 
#         The datasets are given in ‘strike’ and ‘dip’ of the fault plane and 
#         ‘trend’ and ‘plunge’ of the slip-lineation using the convention illustrated 
#         in Fig. 1, i.e. with the direction of strike defined 90&deg; counter-clockwise 
#         to the dip direction. Note that negative plunge represent reverse slip 
#         in the direction opposite the given trend. To use the datasets in the 
#         FSinvert script, copy and paste the data excluding the header into a 
#         .txt file and assign the correct path to the dataset in the script.
#         
# Output: principal stress axes (trend/plunge), shape ratio phi, vorticity axis
#         and magnitude — for both the 9D and 6D formulations of Angelier (1990).

#' 9D Direct Inversion for Fault Slip Including Vorticity 
#' 
#' Direct inversion of stress, strain or strain rate including vorticity using 
#' 9D parameter space using the method by Hansen (2013). It can be applied 
#' regardless whether the dynamic or the kinematic hypothesis is adopted; 
#' it can handle datasets representing two to seven degrees of freedom; and it 
#' is not dependent on the correct assessment of slip sense. 
#'
#' @param x object of class `"Pair"` or `"Fault"` with at least 7 rows.
#' @param flip logical. Flip if you want to have the negative stress tensor, i.e. 
#' sigma 1 and 3 will be flipped.
#' @inheritParams slip_inversion_michael
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
#' par(mfrow = c(1, 2))
#' invisible(lapply(angelier1990, function(x){
#' 
#' res <- slip_inversion_hansen(x, TRUE)
#' phi_val <- round(res$stress_shape$phi, 2)
#' rup_val <-  round(res$misfit$rup, 2)
#' w_val <- round(res$vorticity_mag, 2)
#' 
#' plot(x, col = 'lightgrey')
#' points(res$principal_axes, col = 2:4, pch = 16, cex = 2)
#' text(res$principal_axes, labels = rownames(res$principal_axes), col = 2:4, adj = -.5)
#' points(res$vorticity_axis, col = 5, pch = 17, cex = 2)
#' text(res$vorticity_axis, labels = bquote(omega), col = 5, adj = -.5)
#' 
#' title(sub = bquote(Phi == .(phi_val) ~ "|" ~ bar("RUP") == .(rup_val) * '%'~"|"~omega == .(w_val)))
#' }))
slip_inversion_hansen <- function(x, flip = FALSE, friction = 0.6) {
  tsign <- if(flip) -1 else 1
  
  stopifnot(is.Pair(x))
  normals <- Vec3(Plane(x))
  slips <- if(is.Fault(x)) Ray(x) else Line(x)
  slips <- Vec3(slips)
  
  n <- nrow(normals)
  if (n < 7L)
    warning("Fewer than 7 fault-slip measurements: solution may be ",
            "underdetermined.")
  
  # Row-wise unit normalisation
  # normals <- normals / sqrt(rowSums(normals^2))
  # slips   <- slips   / sqrt(rowSums(slips^2))
  
  # 1. Compute poles to M-planes and f-poles  (Eqs. 14, 19)
  # n_vec = fault plane normal  (n in paper)
  # b_vec = pole to M-plane     (b in paper) = n x v  (Fig. 1)
  # f-pole (bf, Eq. 19): bf = [b1n1, b1n2, b1n3, b2n1, ..., b3n3] / |.|
  
 #  b <- cbind(
 #   normals[, 2] * slips[, 3] - normals[, 3] * slips[, 2],
 #   normals[, 3] * slips[, 1] - normals[, 1] * slips[, 3],
 #   normals[, 1] * slips[, 2] - normals[, 2] * slips[, 1]
 # )
 b <- crossprod(normals, slips)

 f <- cbind(
   b[, 1] * normals[, 1], b[, 1] * normals[, 2], b[, 1] * normals[, 3],
   b[, 2] * normals[, 1], b[, 2] * normals[, 2], b[, 2] * normals[, 3],
   b[, 3] * normals[, 1], b[, 3] * normals[, 2], b[, 3] * normals[, 3]
 )
 f <- f / sqrt(rowSums(f^2))
  
  # 2. Second moment tensor  (Eq. 21)
  M <- crossprod(f)
  
  # 3. Eigen-decomposition; second eigenvector = best estimate of the slip vector
  #    (Eq. 20)
  eM    <- eigen(M, symmetric = TRUE)
  MSort <- order(Re(eM$values))              # ascending
  val_M <- Re(eM$values)[MSort]
  vec_M <- t(Re(eM$vectors))[MSort, ]       # rows = eigenvectors
  #val_M <- Re(eM$values)
  #vec_M <- Re(eM$vectors)
  
  s <- vec_M[2, ]   # second eigenvector in ascending rank
  
  # 4. Build inverted slip tensor; symmetric + antisymmetric decomposition
  #    (Eqs. 8–10; algorithm steps 6–7)
  T_mat      <- matrix(s, 3, 3, byrow = TRUE) * tsign
  Ts         <- (T_mat + t(T_mat)) / 2   # b_T_hat_S
  Ta         <- (T_mat - t(T_mat)) / 2   # b_T_hat_A
  
  # 5. Principal axes and shape ratio from symmetric part  (Eqs. 2–5)
  eTs    <- eigen(Ts, symmetric = TRUE)
  #TsSort <- order(Re(eTs$values), decreasing = TRUE)
  #val_Ts <- Re(eTs$values)[TsSort]
  val_Ts <- eTs$values 
  #vec_Ts <- t(Re(eTs$vectors))[TsSort, ]   # rows = eigenvectors (s1, s2, s3)
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
  ratio          <- Norm / Ts
  ratio[Ts == 0 & Norm == 0] <- 0          # guard 0/0 -> 0
  Ta             <- Ta * ratio
  
  # Axial vector of T_A (Eq. 10)
  w <- as.Vec3(c(Ta[3, 2], Ta[1, 3], Ta[2, 1]))
  w_xyz <- vnorm(w)
  w_mag <- -2 * vlength(w)
  
  # Enforce downward-pointing convention
  #if (nw[3] >= 0) nw <- -nw else wlen <- -wlen
  

  # 7. Convert to [Trend, Plunge]  (Hansen's atan + quadrant logic)
  # vec_Ts <- vec_Ts * ifelse(vec_Ts[, 3] >= 0, -1, 1)
  # Stp <- cbind(
  #   (atan2(vec_Ts[, 1], vec_Ts[, 2]) * 180 / pi) %% 360,
  #   (asin(-vec_Ts[, 3]) * 180 / pi)
  # )
  # Wtp <- .tp(nw)
  principal_axes <- as.Vec3(t(vec_Ts))
  rownames(principal_axes) <- names(val_Ts) <- c('sigma1', 'sigma2', 'sigma3')
  
  
  # 8. 
  
  #p <- tau2stress(Ts)
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

  list(
      stress_tensor = as.ellipsoid(Ts),
      principal_axes = Line(principal_axes),
      principal_vals = val_Ts,
      principal_faults = pfaults,
      
      vorticity_mag = w_mag,
      vorticity_axis = Line(w_xyz),
      
      stress_shape = stress_shape,
      tau_mean = sigma_s_mean,
      stress_components = cbind(shearnorm, tendency),
      misfit = misfit,
      SHmax = SHmax
  )
}

# ---- Internal helpers -------------------------------------------------------

# Convert a downward-pointing unit direction cosine vector to [Trend, Plunge].
# Replicates Hansen's atan() + quadrant-correction block exactly.
# .tp <- function(v) {
#   c((atan2(v[1], v[2]) * 180 / pi),
#     asin(-v[3]) * 180 / pi)
# }
