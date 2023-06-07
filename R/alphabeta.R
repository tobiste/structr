#' Orientation of structures from drill core orientation angles
#'
#' Calculates the orientation of a plane or line from internal core angles (
#' alpha, beta and gamma) of oriented drill cores
#'
#' @param alpha numeric vector. Alpha angle in degrees
#' @param beta numeric vector. Beta angle in degrees
#' @param gamma numeric. Gamma angle in degrees
#' @param azi numeric. Azimuth of drill core axis orientation in degrees (measured clockwise from North).
#' @param inc numeric. Inclination of drill core axis in degrees (Note: **negative values for downward**).
#' @return numeric two-column matrix with Dip direction and dip angle of plane or
#' azimuth and plunge of line (in degrees).
#' @name drillcore
#' @examples
#' azi <- 225
#' inc <- -45
#' plane_from_drillcore(60, 320, azi, inc)
#' plane_from_drillcore(45, 220, azi, inc)
#'
#' # multiple alpha-beta measurements
#' stereoplot()
#' stereo_point(Line(azi, -inc), lab = "CA")
#' plane_from_drillcore(alpha = c(60, 45), beta = c(320, 220), azi, inc) |> 
#' stereo_point(lab = c("A", "B"))
NULL

#' @rdname drillcore
#' @export
plane_from_drillcore <- function(alpha, beta, azi, inc) {
  stopifnot(length(alpha) == length(beta))
  inc <- -inc

  n <- length(alpha)
  if (length(azi) == 1) {
    azi <- rep(azi, n)
    inc <- rep(inc, n)
  }

  C <- line2vec(cbind(azi, inc)) # core axis
  B <- line2vec(cbind((azi + 180) %% 360, 90 - inc)) # bottom of hole line

  P <- E <- int <- matrix(nrow = n, ncol = 3)
  for (i in 1:n) {
    if (beta[i] %% 180 < 90) {
      E[i, ] <- vrotate(B[i, ], C[i, ], -tectonicr::deg2rad(beta[i])) |> vcross(C[i, ])
      int[i, ] <- vcross(E[i, ], C[i, ])
    } else {
      E[i, ] <- vrotate(B[i, ], C[i, ], tectonicr::deg2rad(beta[i]))
      int[i, ] <- vcross(C[i, ], E[i, ])
    }
  }
  P <- vrotate(C, int, -tectonicr::deg2rad(90 - alpha))
  vec2plane(P)
}

#' @rdname drillcore
#' @export
line_from_drillcore <- function(alpha, beta, gamma, azi, inc) {
  stopifnot(length(alpha) == length(beta))
  inc <- -inc

  n <- length(alpha)
  if (length(azi) == 1) {
    azi <- rep(azi, n)
    inc <- rep(inc, n)
  }

  C <- line2vec(cbind(azi, inc)) # core axis
  B <- line2vec(cbind((azi + 180) %% 360, 90 - inc)) # bottom of hole line

  L <- P <- E <- int <- matrix(nrow = n, ncol = 3)
  for (i in 1:n) {
    if (beta[i] %% 180 < 90) {
      E[i, ] <- vrotate(B[i, ], C[i, ], -tectonicr::deg2rad(beta[i])) |> vcross(C[i, ])
      int[i, ] <- vcross(E[i, ], C[i, ])
    } else {
      E[i, ] <- vrotate(B[i, ], C[i, ], tectonicr::deg2rad(beta[i]))
      int[i, ] <- vcross(C[i, ], E[i, ])
    }
  }
  P <- vrotate(C, int, -tectonicr::deg2rad(90 - alpha))
  L <- vrotate(int, P, tectonicr::deg2rad(gamma))
  vec2line(L)
}

#
# azi = 225; inc = -45
# alpha = c(60, 45); beta = c(320, 220)
# gamma = c(10, 20)
# inc <- -inc
#
# n <- length(alpha)
# if(length(azi) == 1){
#   azi <- rep(azi, n)
#   inc <- rep(inc, n)
# }
#
# C <- line2vec(cbind(azi, inc)) # core axis
# B <- line2vec(cbind((azi + 180) %% 360, 90 - inc)) # bottom of hole line
#
#   L <- P <- E <- int <- matrix(nrow = n, ncol = 3)
#   for(i in 1:n){
#     if(beta[i] %% 180 < 90){
#       E[i, ] <- vrotate(B[i, ], C[i, ], -tectonicr::deg2rad(beta[i])) |> vcross(C[i, ])
#       int[i, ] <- vcross(E[i, ], C[i, ])
#     } else {
#       E[i, ] <- vrotate(B[i, ], C[i, ], tectonicr::deg2rad(beta[i]))
#       int[i, ] <- vcross(C[i, ], E[i, ])
#     }
#
#   }
#   P <- vrotate(C , int , -tectonicr::deg2rad(90 - alpha))
#   L <- vrotate(int, P, tectonicr::deg2rad(gamma))
#
#   stereoplot(col = NULL)
#   stereo_point(vec2line(C), lab = "C")
#   stereo_greatcircle(vec2plane(C)[1, ])
#
#   stereo_smallcircle(vec2line(C)[1, ], d = 90-alpha[1], col = 2, lty = 2)
#   stereo_smallcircle(vec2line(C)[2, ], d = 90-alpha[2], col = 3, lty = 2)
#
#   stereo_point(vec2line(B), lab = "B")
#
#   stereo_point(vec2line(E), lab = c("E1", "E2"), col = c(2, 3), pch = 21)
#
#   stereo_point(vec2line(int), lab = c("int1", 'int2'), col = c(2, 3), pch = 22)
#   stereo_greatcircle(vec2plane(int)[1, ], col = c(2), lty = 3)
#   stereo_greatcircle(vec2plane(int)[2, ], col = c(3), lty = 3)
#
#   stereo_point(vec2plane(P), col = c(4, 5), plane = TRUE, lab = c("P1", 'P2'), pch = 23)
#   stereo_greatcircle(vec2plane(P)[1, ], col = c(4), lty = 4)
#   stereo_greatcircle(vec2plane(P)[2, ], col = c(5), lty = 4)
#
#   stereo_point(vec2line(L), col = c(4, 5), plane = F, lab = c("L1", 'L2'), pch = 24)
#
