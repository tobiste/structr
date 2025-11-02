#' Direction of maximum horizontal stress from the stress tensor
#'
#' Calculates the direction of maximum horizontal stress using only the
#' directions of the principal stress and \eqn{R = \frac{S1 - S2}{S1 - S3}}.
#' This function Equations 11 and 10 from Lund and Townend (2007).
#'
#' @param S1,S2,S3 The principal stress orientations.
#' The variables hold the coordinates in the North, East and Down geographical
#' coordinate system, e.g. `S1 = c(s1N,s1E,s1D)`. Given as object of class `"Vec3"` or `"Line"`
#' @param R numeric. Relative magnitude of `S2` with respect to `S1` and `S3`:
#' \eqn{R = \frac{S1 - S2}{S1 - S3}}. Values ranging from 0 to 1, with 0 being
#' `S1==S2` and 1 being `S2==S3`. Equivalent to the stress shape ratio of Gephart & Forsyth (1984).
#' @param tol Tolerance of comparison.
#' @param ortho.tol tolerance angle (in degree) for orthogonality check of the
#' three principal stress vectors.
#'
#' @return The direction of SH from North as numeric angle in degrees (radians if all principal stress axes were given as `"Vec3"` objects).
#' @export
#'
#' @seealso [SH_from_tensor()] when full knowledge of stress tensor; [slip_inversion()] for stress inversion of fault slip data.
#'
#' @references Lund and Townend, (2007). Calculating horizontal stress
#' orientations with full or partial knowledge of the tectonic stress tensor,
#' Geophys. J. Int., doi:\doi{10.1111/j.1365-246X.2007.03468.x}.
#'
#' @examples
#' # first example from https://www.snsn.se/SH/SHcode/benchmark.out
#' S1 <- Line(250.89, 70.07)
#' S3 <- Line(103.01, 17.07)
#' S2 <- crossprod(S3, S1)
#' SH(S1, S2, S3, R = 1) #  70.89
#'
#' R <- seq(0, 1, .05)
#' cbind(R, SH = SH(S1, S2, S3, R = R))
SH <- function(S1, S2, S3, R, tol = .Machine$double.eps^0.5, ortho.tol = 0.005) {
  # Convert to unit vectors
  all_rad <- all(is.Vec3(S1), is.Vec3(S2), is.Vec3(S3))

  S1 <- Vec3(S1) |> vnorm()
  S2 <- Vec3(S2) |> vnorm()
  S3 <- Vec3(S3) |> vnorm()

  n <- nrow(S1)
  stopifnot(n == nrow(S2), n == nrow(S3))

  # Check orthogonality
  ortho.tol <- deg2rad(ortho.tol)
  if (dotprod(S1, S2) > ortho.tol) stop("S1 and S2 are not orthogonal.")
  if (dotprod(S1, S3) > ortho.tol) stop("S1 and S3 are not orthogonal.")
  if (dotprod(S2, S3) > ortho.tol) stop("S2 and S3 are not orthogonal.")

  # Check R validity
  if (length(R) == 1 & n > 1) {
    R <- rep(R, n)
  }

  S1 <- unclass(S1)
  S2 <- unclass(S2)
  S3 <- unclass(S3)

  if (!is.numeric(R) || any(R < 0) || any(R > 1)) stop("SH(): R must be between 0 and 1.")

  # Calculate the denominator and numerator in Eq. 11
  X <- S1[, 1] * S1[, 1] - S1[, 2] * S1[, 2] + (1 - R) * (S2[, 1] * S2[, 1] - S2[, 2] * S2[, 2])
  Y <- 2 * (S1[, 1] * S1[, 2] + (1 - R) * (S2[, 1] * S2[, 2]))
  # print 'SH: X %f Y %f' % (X,Y)

  # Is the denominator (X here) from Eq. 11 in Lund and Townend (2007) zero?

  alpha <- vapply(seq_along(R), function(i) {
    if (abs(X[i]) < tol) {
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

      if (abs(S2[i, 1] * S2[i, 2]) < tol) {
        # s2Ns2E = 0
        # print 'SH: s2Ns2E = 0'

        # print 'SH: s1Ns1E = s2Ns2E = 0 => SH undefined'
        if (abs(S1[i, 1] * S1[i, 2]) < tol) NA_real_ else pi / 4
      } else {
        # print 'SH: R fits => SH undefined'
        if (abs(R[i] - (1 + S1[i, 1] * S1[i, 2] / S2[i, 1] * S2[i, 2])) < tol) NA_real_ else pi / 4
      }
    } else {
      # The denominator is non-zero
      # print 'SH: X != 0 Denominator non-zero'
      # alpha <- atan(Y / X) / 2
      atan2(Y[i], X[i]) / 2
    }
  }, FUN.VALUE = numeric(n))


  # Have we found a minimum or maximum? Use 2nd derivative to find out.
  # A negative 2nd derivative indicates a maximum, which is what we want.
  # print X,Y,alpha
  dev2 <- -2 * X * cos(2 * alpha) - 2 * Y * sin(2 * alpha)
  # print 'SH: dev2 %f' % (dev2)

  # We found a minimum. Add 90 degrees to get the maximum.
  # if (dev2 > 0) alpha <- alpha + pi / 2
  alpha <- ifelse(dev2 > 0, alpha + pi / 2, alpha)

  # The resulting direction of SH is given as [0,180[ degrees.

  # if (alpha < 0) {
  #   alpha <- alpha + pi
  # }
  # if (alpha > pi | abs(alpha - pi) < tol) {
  #   alpha <- alpha - pi
  # }
  alpha <- alpha %% pi

  # Return angle in radians only if all principles axes were given as Cartesian vectors
  if (!all_rad) alpha <- rad2deg(alpha)


  return(unname(alpha))
}

#' Direction of maximum horizontal stress from the stress tensor (full knowledge)
#'
#' @param S 3x3 matrix where the columns are the principal stress axes and the rows are the coordinates.
#'
#' @return numeric angle in degrees. The direction of SH from North.
#' @export
#'
#' @references Lund and Townend, (2007). Calculating horizontal stress
#' orientations with full or partial knowledge of the tectonic stress tensor,
#' Geophys. J. Int., doi:\doi{10.1111/j.1365-246X.2007.03468.x}.
#' 
#' @seealso [SH()] when only principal axes and their relative magnitudes are known; 
#' [slip_inversion()] for stress inversion of fault slip data.
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

  # rownames(S) <- c('x', 'y', 'z')

  SH(
    as.Vec3(S[, 1]),
    as.Vec3(S[, 2]),
    as.Vec3(S[, 3]),
    R = (S1 - S2) / (S1 - S3)
  ) |>
    rad2deg()

  # Y <- 2 * (S1 * S[1, 1] * S[2, 1] + S2 * S[1, 2] * S[2, 2] + S3 * S[1, 3] * S[2, 3])
  # X <- S1 * (S[1, 1]^2 - S[2, 1]^2) + S2 * (S[1, 2]^2 - S[2, 2]^2) + S3 * (S[1, 3]^2 - S[2, 3]^2)
  #
  # # alpha0 <- rad2deg(atan(Y / X) / 2)
  # alpha <- rad2deg(atan2(Y, X) / 2) %% 180
  # return(alpha)
}
