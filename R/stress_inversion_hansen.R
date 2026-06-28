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
#' res <- slip_inversion_hansen(x, TRUE)
#' plot(x, col = 'lightgrey')
#' title(sub = paste0("R = ", round(res$phi, 2)))
#' points(res$principal_axes, col = 2:4, pch = 16, cex = 2)
#' text(res$principal_axes, labels = 1:3, col = 2:4, adj = -1)
#' points(res$vorticity_axis, col = 5, pch = 17, cex = 2)
#' }))
slip_inversion_hansen <- function(x, flip = FALSE) {
  tsign <- if(flip) -1 else 1
  
  stopifnot(is.Pair(x))
  normals <- unclass(Vec3(Plane(x)))
  slips <- if(is.Fault(x)) Ray(x) else Line(x)
  slips <- unclass(Vec3(slips))
  
  n <- nrow(normals)
  if (n < 7L)
    warning("Fewer than 7 fault-slip measurements: solution may be ",
            "underdetermined (Hansen 2013, §4).")
  
  # Row-wise unit normalisation
  # normals <- normals / sqrt(rowSums(normals^2))
  # slips   <- slips   / sqrt(rowSums(slips^2))
  
  # --------------------------------------------------------------------------
  # 1. Compute poles to M-planes and f-poles  (Eqs. 14, 19)
  # --------------------------------------------------------------------------
  # n_vec = fault plane normal  (n in paper)
  # b_vec = pole to M-plane     (b in paper) = n x v  (Fig. 1)
  # f-pole (bf, Eq. 19): bf = [b1n1, b1n2, b1n3, b2n1, ..., b3n3] / |.|
  
  f <- matrix(0, n, 9)
  
  for (i in seq_len(n)) {
    p <- normals[i, ]                      # n in paper
    m <- .cross3(p, slips[i, ])            # b in paper = n x v
    
    fi <- c(m[1]*p[1], m[1]*p[2], m[1]*p[3],
            m[2]*p[1], m[2]*p[2], m[2]*p[3],
            m[3]*p[1], m[3]*p[2], m[3]*p[3])
    f[i, ] <- fi / sqrt(sum(fi * fi))
  }
  
  # --------------------------------------------------------------------------
  # 2. Second moment tensor  (Eq. 21)
  # --------------------------------------------------------------------------
  M <- matrix(0, 9, 9)
  for (i in seq_len(n)) M <- M + f[i, ] %o% f[i, ]
  
  # --------------------------------------------------------------------------
  # 3. Eigen-decomposition; second-lowest eigenvector = best estimate of bs
  #    (Eq. 20)
  # --------------------------------------------------------------------------
  eM    <- eigen(M)
  MSort <- order(Re(eM$values))              # ascending
  val_M <- Re(eM$values)[MSort]
  vec_M <- t(Re(eM$vectors))[MSort, ]       # rows = eigenvectors
  
  s <- vec_M[2, ]   # second eigenvector in ascending rank
  
  # --------------------------------------------------------------------------
  # 4. Build inverted slip tensor; symmetric + antisymmetric decomposition
  #    (Eqs. 8–10; algorithm steps 6–7)
  # --------------------------------------------------------------------------
  T_mat      <- matrix(s, 3, 3, byrow = TRUE) * tsign
  Ts         <- (T_mat + t(T_mat)) / 2   # b_T_hat_S
  Ta         <- (T_mat - t(T_mat)) / 2   # b_T_hat_A
  
  # --------------------------------------------------------------------------
  # 5. Principal axes and shape ratio from symmetric part  (Eqs. 2–5)
  # --------------------------------------------------------------------------
  eTs    <- eigen(Ts, symmetric = TRUE)
  #TsSort <- order(Re(eTs$values), decreasing = TRUE)
  #val_Ts <- Re(eTs$values)[TsSort]
  val_Ts <- eTs$values 
  #vec_Ts <- t(Re(eTs$vectors))[TsSort, ]   # rows = eigenvectors (s1, s2, s3)
  vec_Ts <- t(eTs$vectors)  
  
  nval_Ts <- (val_Ts - val_Ts[3]) / (val_Ts[1] - val_Ts[3])
  phi     <- nval_Ts[2]
  
  # Reconstruct reduced symmetric tensor T_S = V D V^T  (Eq. 5)
  # V has eigenvectors as columns = t(vec_Ts)
  Norm <- t(vec_Ts) %*% diag(c(nval_Ts[1], nval_Ts[2], 0)) %*% vec_Ts
  
  # --------------------------------------------------------------------------
  # 6. Vorticity axis and magnitude  (Eqs. 12–13, 22)
  # --------------------------------------------------------------------------
  # Normalise antisymmetric part: b_T_A = b_T_hat_A * (T_S / b_T_hat_S)
  ratio          <- Norm / Ts
  ratio[Ts == 0 & Norm == 0] <- 0          # guard 0/0 -> 0
  Ta             <- Ta * ratio
  
  # Axial vector of b_T_A (Eq. 10 sign layout; matches Python index convention)
  w    <- c(Ta[3, 2], Ta[1, 3], Ta[2, 1])
  wlen <- sqrt(sum(w * w))
  nw   <- w / wlen
  
  # Enforce downward-pointing convention
  if (nw[3] >= 0) nw <- -nw else wlen <- -wlen
  
  # --------------------------------------------------------------------------
  # 7. Convert to [Trend, Plunge]  (Hansen's atan + quadrant logic)
  # --------------------------------------------------------------------------
  Stp <- matrix(0, 3, 2)
  for (i in 1:3) {
    if (vec_Ts[i, 3] >= 0) vec_Ts[i, ] <- -vec_Ts[i, ]
    Stp[i, ] <- .tp(vec_Ts[i, ])
  }
  Wtp <- .tp(nw)
  
  # --------------------------------------------------------------------------
  # 8. Print results
  # --------------------------------------------------------------------------
  
    list(
      stress_tensor = as.ellipsoid(Ts),
      principal_axes = as.Ray(Stp),
      principal_vals = val_Ts,
      phi = phi,
      #full = T_mat,
      #assymmetric = as.ellipsoid(Ta),
      vorticity_mag = wlen * 2 * tsign,
      vorticity_axis = as.Ray(Wtp)
  )

}

# ---- Internal helpers -------------------------------------------------------

.cross3 <- function(a, b) {
  c(a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1])
}

# Convert a downward-pointing unit direction cosine vector to [Trend, Plunge].
# Replicates Hansen's atan() + quadrant-correction block exactly.
.tp <- function(v) {
  r2d <- 180 / pi
  trend  <- (atan(v[1] / v[2]) * r2d)
  if      (v[1] >= 0 && v[2] <  0) trend <- trend + 180
  else if (v[1] <  0 && v[2] <  0) trend <- trend + 180
  else if (v[1] <  0 && v[2] >= 0) trend <- trend + 360
  plunge <- (asin(-v[3]) * r2d)
  c(trend, plunge)
}

