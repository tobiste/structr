#' Orientation of structures from drill core orientation angles
#'
#' Calculates the orientation of a plane or line from alpha, beta and gamma
#' angles in oriented drill cores
#'
#' @param alpha numeric vector. Alpha angle in degrees
#' @param beta numeric vector. Beta angle in degrees
#' @param gamma numeric. Gamma angle in degrees
#' @param azi numeric. Azimuth of drill core axis orientation in degrees (measured from North).
#' @param inc numeric. Inclination of drill core axis in degrees (negative values for downward!).
#' @return numeric two-column matrix with Dip direction and dip angle of plane or
#' azimuth and plunge of line (in degrees).
#' @name drillcore
#' @examples
#' drill_plane(45, 220, 225, -45)
#' drill_plane(60, 320, 225, -45)
#'
#' # multiple alpha-beta measurements
#' drill_plane(alpha = c(45, 60), beta = c(220, 320), azi = 225, inc = -45)
NULL

#' @rdname drillcore
#' @export
drill_plane <- function(alpha, beta, azi, inc) {
  stopifnot(length(alpha) == length(beta), length(azi) == 1, length(inc) == 1)
  n <- length(alpha)
  azi <- rep(azi, n)
  inc <- rep(-inc, n)

  C <- line2vec(cbind(azi, inc))
  B <- line2vec(cbind(azi + 180, 90 - inc))

  E <- vrotate(B, C, beta * DEG2RAD())
  I <- vcross(C, E)
  P <- vrotate(C, I, pi - alpha * DEG2RAD())
  vec2plane(P)
}

#' @rdname drillcore
#' @export
drill_line <- function(alpha, beta, gamma, azi, inc) {
  inc <- -inc
  C <- line2vec(cbind(azi, inc))
  B <- line2vec(cbind(azi + 180, 90 - inc))

  E <- vrotate(B, C, beta * DEG2RAD())
  I <- vcross(C, E)
  P <- vrotate(C, I, pi - alpha * DEG2RAD())
  L <- vrotate(I, P, pi + gamma * DEG2RAD())
  vec2line(L)
}

# alpha = 45; beta = 220; gamma = 20
# azi = 225; inc = -45
#
# inc <- -inc
# CA <- lin2vec(azi, inc)
# B <- lin2vec(azi + 180, 90-inc)
#
# E <- vrotate(B, CA, tectonicr::deg2rad(beta))
# P <- vrotate(CA, vcross(CA, E), pi - tectonicr::deg2rad(alpha))
#
# I <- vcross(E, P)
# L <- vrotate(I, P, pi + tectonicr::deg2rad(gamma))
#
# CA.lin = vec2lin(CA)
# B.lin = vec2lin(B)
# E.lin = vec2lin(E)
# P.fol <- vec2fol(P)
# I.lin = vec2lin(I)
# L.lin = vec2lin(L)
#
# I2.lin <- vcross(CA, P) |> vec2lin()
# I3.lin <- vcross(CA, E) |> vec2lin()
#
#
# RFOC::net(col = "grey90")
# stereo_plane(data.frame(P_dipdir = I.lin[1]+180, P_dip = 90-I.lin[2]), col = 5)
# stereo_plane(data.frame(P_dipdir = P.fol[1], P_dip = P.fol[2]), col = 4)
#
# stereo_line(data.frame(L_azimuth = CA.lin[1], L_plunge = CA.lin[2]), col = 1, lab = "CA")
# stereo_plane(data.frame(P_dipdir = 180+CA.lin[1], P_dip = CA.lin[2]))
# stereo_line(data.frame(L_azimuth = B.lin[1], L_plunge = B.lin[2]), col = 2, lab = "B")
# stereo_line(data.frame(L_azimuth = E.lin[1], L_plunge = E.lin[2]), col = 3, lab = "E")
# stereo_pole(data.frame(P_dipdir = P.fol[1], P_dip = P.fol[2]), col = 4, lab = "P")
# stereo_line(data.frame(L_azimuth = I.lin[1], L_plunge = I.lin[2]), col = 5, lab = "I")
# stereo_line(data.frame(L_azimuth = L.lin[1], L_plunge = L.lin[2]), col = 6, lab = "L")
