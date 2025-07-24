#' Direction of maximum horizontal stress from the stress tensor
#'
#' Calculates the direction of maximum horizontal stress using only the
#' directions of the principal stress and \eqn{R = \frac{S1 - S2}{S1 - S3}}.
#' This function Equations 11 and 10 from Lund and Townend (2007).
#'
#' @param S1,S2,S3 numeric three-column vectors. The principal stress orientations.
#' The variables hold the coordinates in the North, East and Down geographical
#' coordinate system, e.g. `S1 = c(s1N,s1E,s1D)`
#' @param R numeric. Relative magnitude of `S2` with respect to `S1` and `S3`:
#' \eqn{R = \frac{S1 - S2}{S1 - S3}}. Values ranging from 0 to 1, with 0 being
#' `S1==S2` and 1 being `S2==S3`.
#' @param tol Tolerance of comparison.
#' @param ortho.tol tolerance angle (in radians) for orthogonality check of the
#' three principal stress vectors.
#'
#' @return numeric angle in degrees. The direction of SH from North.
#' @export
#'
#' @references Lund and Townend, (2007). Calculating horizontal stress
#' orientations with full or partial knowledge of the tectonic stress tensor,
#' Geophys. J. Int., doi:\doi{10.1111/j.1365-246X.2007.03468.x}.
#'
#' @examples
#' # first example from https://www.snsn.se/SH/SHcode/benchmark.out
#' S1 <- Line(250.89, 70.07)
#' S3 <- Line(103.01, 17.07)
#' S2 <- vcross(S3, S1)
#' SH(S1, S2, S3, R = 1) #  70.89 
#'
#' R <- seq(0, 1, .05)
#' cbind(R, SH = sapply(R, function(x) {
#'   SH(S1, S2, S3, R = x)
#' }))
SH <- function(S1, S2, S3, R, tol = .Machine$double.eps^0.5, ortho.tol = 1e-4) {
  # Convert to unit vectors
  S1 <- if (is.spherical(S1)) line2vec(S1) else vnorm(vec2mat(S1))
  S2 <- if (is.spherical(S2)) line2vec(S2) else vnorm(vec2mat(S2))
  S3 <- if (is.spherical(S3)) line2vec(S3) else vnorm(vec2mat(S3))
  
  # Check orthogonality
  if (vdot(S1, S2) > ortho.tol) stop("S1 and S2 are not orthogonal.")
  if (vdot(S1, S3) > ortho.tol) stop("S1 and S3 are not orthogonal.")
  if (vdot(S2, S3) > ortho.tol) stop("S2 and S3 are not orthogonal.")
  
  # Check R validity
  if (!is.numeric(R) || R < 0 || R > 1) stop("SH(): R must be between 0 and 1.")
  
  # Calculate the denominator and numerator in Eq. 11
  X <- S1[1] * S1[1] - S1[2] * S1[2] + (1 - R) * (S2[1] * S2[1] - S2[2] * S2[2])
  Y <- 2 * (S1[1] * S1[2] + (1 - R) * (S2[1] * S2[2]))
  # print 'SH: X %f Y %f' % (X,Y)

  # Is the denominator (X here) from Eq. 11 in Lund and Townend (2007) zero?
  if (abs(X) < tol) {
    # print 'SH: X = 0 R = %.2f' % R
    # If so, the first term in Eq. 10 is zero and we are left with the
    # second term, which must equal zero for a stationary point.
    # The second term is zero either if
    # s1Ns1E + (1-R)*s2Ns2E = 0   (A)
    # or if
    # cos(2*alpha) = 0            (B)
    # If (A) holds, the direction of SH is undefined since Eq. 10 is zero
    # irrespective of the value of alpha. We therefore check for (A) first.
    # If (A) holds, R = 1 + s1Ns1E/s2Ns2E unless s2Ns2E = 0, in which case
    # s1Ns1E also has to be zero for (A) to hold.

    if (abs(S2[1] * S2[2]) < tol) {
      # s2Ns2E = 0
      # print 'SH: s2Ns2E = 0'
      
      # print 'SH: s1Ns1E = s2Ns2E = 0 => SH undefined'
      if (abs(S1[1] * S1[2]) < tol) alpha <- NA_real_ else alpha <- pi / 4
    } else {
      # print 'SH: R fits => SH undefined'
      if (abs(R - (1 + S1[1] * S1[2] / S2[1] * S2[2])) < tol) alpha <- NA_real_ else alpha <- pi / 4
    }
  } else {
    # The denominator is non-zero
    # print 'SH: X != 0 Denominator non-zero'
    # alpha <- atan(Y / X) / 2
    alpha <- atan2(Y, X) / 2
  }

  # Have we found a minimum or maximum? Use 2nd derivative to find out.
  # A negative 2nd derivative indicates a maximum, which is what we want.
  # print X,Y,alpha
  dev2 <- -2 * X * cos(2 * alpha) - 2 * Y * sin(2 * alpha)
  # print 'SH: dev2 %f' % (dev2)
  
  # We found a minimum. Add 90 degrees to get the maximum.
  if (dev2 > 0) alpha <- alpha + pi / 2
  
  # The resulting direction of SH is given as [0,180[ degrees.
  
  # if (alpha < 0) {
  #   alpha <- alpha + pi
  # }
  # if (alpha > pi | abs(alpha - pi) < tol) {
  #   alpha <- alpha - pi
  # }
  alpha <- alpha %% pi
  
  return(rad2deg(alpha))
}

#' Direction of maximum horizontal stress from the stress tensor (full knowledge)
#'
#' @param S 3x3 matrix
#'
#' @return numeric angle in degrees. The direction of SH from North.
#' @export
#'
#' @references Lund and Townend, (2007). Calculating horizontal stress
#' orientations with full or partial knowledge of the tectonic stress tensor,
#' Geophys. J. Int., doi:\doi{10.1111/j.1365-246X.2007.03468.x}.
#'
#' @examples
#' S1 <- c(-0.11163471165673433, -0.32215903156083286, 0.94008044843891103)
#' S2 <- c(0.97015303285724142, 0.16959307992318776, 0.17332420511879815)
#' S3 <- c(-0.21526909669344957, 0.93137089584437482, 0.29361108695489119)
#'
#' S <- cbind(S1, S2, S3) * cbind(rep(3, 3), rep(2, 3), rep(1, 3)) / 3
#' SH_from_tensor(S)
SH_from_tensor <- function(S) {
  SM <- eigen(S %*% t(S), TRUE)
  S1 <- SM$values[1]
  S2 <- SM$values[2]
  S3 <- SM$values[3]

  SH(S[, 1], S[, 2], S[, 3], R = (S1 - S2) / (S1 - S3))

  # Y <- 2 * (S1 * S[1, 1] * S[2, 1] + S2 * S[1, 2] * S[2, 2] + S3 * S[1, 3] * S[2, 3])
  # X <- S1 * (S[1, 1]^2 - S[2, 1]^2) + S2 * (S[1, 2]^2 - S[2, 2]^2) + S3 * (S[1, 3]^2 - S[2, 3]^2)
  #
  # # alpha0 <- rad2deg(atan(Y / X) / 2)
  # alpha <- rad2deg(atan2(Y, X) / 2) %% 180
  # return(alpha)
}
