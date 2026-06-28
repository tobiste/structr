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

#' 9D Direct Inversion of Fault Slip 
#' 
#' Direct inversion of stress, strain or strain rate including vorticity using 
#' 9D parameter space using the method by Hansen (2013)
#'
#' @param x object of class `"Pair"` or `"Fault"`
#' @param flip logical. Flip if you want to have the negative stress tensor, i.e. sigma 1 and 3 will be flipped.
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
#' res <- slip_inversion_hansen(x)
#' plot(x)
#' title(sub = paste0("R = ", round(res$nine$phi, 2)))
#' points(res$nine$principal_axes, col = 2:4, pch = 16, cex = 2)
#' text(res$nine$principal_axes, labels = rownames(res$nine$principal_axes), col = 2:4, adj = -0.5)
#' points(res$six$principal_axes, col = 2:4, pch = 15, cex = 1)
#' points(res$nine$vorticity_axis, col = 5, pch = 17, cex = 2)
#' }))
slip_inversion_hansen <- function(x, flip = FALSE) {
  tsign <- if(flip) -1 else 1
  
  
  stopifnot(is.Pair(x))
  
  normals <- Plane(x)
  slips <- if(is.Fault(x)) Ray(x) else Line(x)
  
  # ---- 1. Read data ---------------------------------------------------------
  
  # inn <- unclass(x)
  n   <- nrow(x)
  
  # ---- 2. Allocate arrays ---------------------------------------------------
  ddr    <- matrix(0, n, 2)   # dip direction, dip
  ldr    <- matrix(0, n, 2)   # lineation trend, plunge
  dcos_p <- matrix(0, n, 3)   # direction cosines of pole to plane
  dcos_l <- matrix(0, n, 3)   # direction cosines of lineation
  dcos_m <- matrix(0, n, 3)   # direction cosines of pole to M-plane
  f      <- matrix(0, n, 9)   # 9D f-poles
  f6     <- matrix(0, n, 6)   # 6D f-poles
  
  # ---- 3. Direction cosines and f-poles -------------------------------------
  deg2rad <- pi / 180
  
  # Dip direction and dip of fault plane
  # ddr[i, 1] <- (inn[i, 1] + 90) %% 360
  #ddr[i, 1] <- inn[i, 1]
  #ddr[i, 2] <- inn[i, 2]
  
  # Lineation
  #ldr[i, 1] <- inn[i, 3]
  #ldr[i, 2] <- inn[i, 4]
  
  dd <- normals[, 1] * deg2rad
  dp <- normals[, 2] * deg2rad
  lt <- slips[, 1] * deg2rad
  lp <- slips[, 2] * deg2rad
  
  for (i in seq_len(n)) {
    # Pole to plane (downward-pointing normal)
    dcos_p[i, ] <- c(
      -sin(dp[i]) * sin(dd[i]),
      -sin(dp[i]) * cos(dd[i]),
      -cos(dp[i])
    )
    
    # Lineation vector
    dcos_l[i, ] <- c(
      cos(lp[i]) * sin(lt[i]),
      cos(lp[i]) * cos(lt[i]),
      -sin(lp[i])
    )
    
    # Pole to M-plane (cross product of pole and lineation)
    dcos_m[i, ] <- .cross3(dcos_p[i, ], dcos_l[i, ])
    
    p <- dcos_p[i, ]
    m <- dcos_m[i, ]
    
    # 9D f-pole (outer product, vectorised, then normalised)
    fi <- c(
      m[1]*p[1], m[1]*p[2], m[1]*p[3],
      m[2]*p[1], m[2]*p[2], m[2]*p[3],
      m[3]*p[1], m[3]*p[2], m[3]*p[3]
    )
    f[i, ] <- fi / sqrt(sum(fi^2))
    
    # 6D f-pole (Voigt-like symmetric reduction, then normalised)
    f6i <- c(
      m[1]*p[1],
      m[2]*p[2],
      m[3]*p[3],
      m[1]*p[2] + m[2]*p[1],
      m[2]*p[3] + m[3]*p[2],
      m[1]*p[3] + m[3]*p[1]
    )
    f6[i, ] <- f6i / sqrt(sum(f6i^2))
  }
  
  # ---- 4. Second moment tensors ---------------------------------------------
  M  <- crossprod(f)    # t(f) %*% f  — same as summing outer products
  M6 <- crossprod(f6)
  
  # ---- 5. Eigen-decomposition -----------------------------------------------
  eM  <- eigen(M,  symmetric = TRUE)   # returns values in *decreasing* order
  eM6 <- eigen(M6, symmetric = TRUE)
  
  # Reorder ascending (Python argsort default), then take index [2] (middle)
  ord  <- order(eM$values)
  ord6 <- order(eM6$values)
  
  s  <- eM$vectors[, ord[2]]    # eigenvector for median eigenvalue
  s6 <- eM6$vectors[, ord6[2]]
  
  # ---- 6. Build stress tensors ----------------------------------------------
  T_mat <- matrix(s, 3, 3, byrow = TRUE) * tsign
  
  Ts <- (T_mat + t(T_mat)) / 2   # symmetric part
  Ta <- (T_mat - t(T_mat)) / 2   # antisymmetric part
  
  T6_mat <- matrix(0, 3, 3)
  T6_mat[1, 1] <- s6[1]; T6_mat[2, 2] <- s6[2]; T6_mat[3, 3] <- s6[3]
  T6_mat[1, 2] <- T6_mat[2, 1] <- s6[4]
  T6_mat[2, 3] <- T6_mat[3, 2] <- s6[5]
  T6_mat[1, 3] <- T6_mat[3, 1] <- s6[6]
  T6_mat <- T6_mat * tsign
  
  # ---- 7. Principal axes and shape ratio ------------------------------------
  eTs <- eigen(Ts, symmetric = TRUE)
  #oTs <- order(eTs$values, decreasing = TRUE)
  val_Ts <- eTs$values
  vec_Ts <- t(eTs$vectors)  # rows = eigenvectors
  
  nval_Ts <- (val_Ts - val_Ts[3]) / (val_Ts[1] - val_Ts[3])
  phi     <- nval_Ts[2]
  
  # Normalised tensor (for vorticity scaling)
  Norm <- vec_Ts %*% diag(c(nval_Ts[1], nval_Ts[2], 0)) %*% t(vec_Ts)
  
  eT6 <- eigen(T6_mat, symmetric = TRUE)
  # oT6 <- order(eT6$values, decreasing = TRUE)
  val_T6 <- eT6$values
  vec_T6 <- t(eT6$vectors)
  
  nval_T6 <- (val_T6 - val_T6[3]) / (val_T6[1] - val_T6[3])
  phi6    <- nval_T6[2]
  
  # ---- 8. Vorticity axis ----------------------------------------------------
  Ta_norm <- Ta * (Norm / Ts)    # element-wise scale (mirrors Python)
  w    <- c(Ta_norm[3, 2], Ta_norm[1, 3], Ta_norm[2, 1])
  wlen <- sqrt(sum(w^2))
  nw   <- w / wlen
  
  # Enforce downward-pointing convention
  if (nw[3] >= 0) {
    nw <- -nw
  } else {
    wlen <- -wlen
  }
  
  # ---- 9. Convert eigenvectors to trend / plunge ----------------------------
  Stp  <- .vec_to_tp(vec_Ts) 
  S6tp <- .vec_to_tp(vec_T6) 
  Wtp  <- .dir_to_tp(nw) 
  rownames(Stp) <- rownames(S6tp) <- c("sigma1", "sigma2", "sigma3")
  
  # ---- 10. Print results ----------------------------------------------------
  # cat(".......................\n\n")
  # cat("9d-space:\n")
  # cat("phi =", round(phi, 2), "\n")
  # cat("s1  =", Stp[1, ], "\n")
  # cat("s2  =", Stp[2, ], "\n")
  # cat("s3  =", Stp[3, ], "\n\n")
  # cat("|vorticity| =", round(wlen * 2, 2), "\n")
  # cat("vorticity   =", Wtp, "\n\n")
  # cat("Eigenvalues:\n"); print(eM$values[ord])
  # cat("\n6d-space:\n")
  # cat("phi =", round(phi6, 2), "\n")
  # cat("s1  =", S6tp[1, ], "\n")
  # cat("s2  =", S6tp[2, ], "\n")
  # cat("s3  =", S6tp[3, ], "\n\n")
  # cat("Eigenvalues:\n"); print(eM6$values[ord6])
  
  # invisible(list(
  #   phi = phi, phi6 = phi6,
  #   s1 = Stp[1, ], s2 = Stp[2, ], s3 = Stp[3, ],
  #   s1_6d = S6tp[1, ], s2_6d = S6tp[2, ], s3_6d = S6tp[3, ],
  #   vorticity_mag = round(wlen * 2, 2),
  #   vorticity_axis = Wtp,
  #   eigenvalues_9d = eM$values[ord],
  #   eigenvalues_6d = eM6$values[ord6]
  # ))
  
  # invisible(
    list(
    nine = list(
      stress_tensor = Ts,
      principal_axes = as.Line(Stp),
      phi = phi,
      principal_vals = val_Ts,
      vorticity_mag = wlen * 2,
      vorticity_axis = as.Line(Wtp)
    ),
    six = list(
      stress_tensor = T6_mat,
      principal_axes = as.Line(S6tp),
      phi = phi6,
      principal_vals = val_T6
    )
  )
  # )
}

# ---- Internal helpers -------------------------------------------------------

# 3-vector cross product
.cross3 <- function(a, b) {
  c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
}

# Convert a direction cosine vector to [trend, plunge] in degrees.
# Enforces downward-pointing (z >= 0 flipped) before conversion.
.dir_to_tp <- function(v) {
  if (v[3] >= 0) v <- -v
  trend <- (atan2(v[1], v[2]) * 180 / pi) %% 360
  plunge <- (asin(-v[3]) * 180 / pi)
  c(trend, plunge)
}

# Apply .dir_to_tp row-wise to a matrix of eigenvectors (one per row).
.vec_to_tp <- function(vecs) {
  tp <- matrix(0, nrow(vecs), 2)
  for (i in seq_len(nrow(vecs))) {
    tp[i, ] <- .dir_to_tp(vecs[i, ])
  }
  tp
}

