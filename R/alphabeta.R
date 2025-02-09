#' Orientation of structures from drill core orientation angles
#'
#' Calculates the orientation of a plane or line from internal core angles (
#' alpha, beta, and gamma) of oriented drill cores
#'
#' @param azi numeric. Azimuth of drill core axis orientation in degrees
#' (measured clockwise from North).
#' @param inc numeric. Inclination of drill core axis in degrees
#' (Note: **negative values for downward**).
#' @param alpha numeric vector. Alpha angle in degrees
#' @param beta numeric vector. Beta angle in degrees
#' @param gamma numeric. (optional). Gamma angle in degrees
#' 
#' @export
#' 
#' @return object of class `"plane"`. If gamma is specified, `"line"` object is
#' returned.
#' 
#' @examples
#' azi <- 225
#' inc <- -45
#' drillcore_orientation(azi, inc, 60, 320)
#' drillcore_orientation(azi, inc, 45, 220)
#'
#' # multiple alpha-beta measurements
#' stereoplot()
#' stereo_point(Line(azi, -inc), lab = "CA")
#' drillcore_orientation(azi, inc, alpha = c(60, 45), beta = c(320, 220)) |>
#'   stereo_point(lab = c("A", "B"))
drillcore_orientation <- function(azi, inc, alpha, beta, gamma = NULL) {
  stopifnot(length(alpha) == length(beta))
  inc <- -inc

  n <- length(alpha)
  if (length(azi) == 1) {
    azi <- rep(azi, n)
    inc <- rep(inc, n)
  }

  C <- line2vec(cbind(azi, inc)) # core axis
  B <- line2vec(cbind((azi + 180) %% 360, 90 - inc)) # bottom of hole line

  E <- int <- matrix(nrow = n, ncol = 3)
  for (i in 1:n) {
    if (cosd(beta[i]) <= 0) {
      E[i, ] <- vrotate(B[i, ], C[i, ], -deg2rad(beta[i])) |> vcross(C[i, ])
      int[i, ] <- vcross(E[i, ], C[i, ])
    } else {
      E[i, ] <- vrotate(B[i, ], C[i, ], deg2rad(beta[i]))
      int[i, ] <- vcross(C[i, ], E[i, ])
    }
  }
  P <- vrotate(C, int, -deg2rad(90 - alpha))

  if (is.null(gamma)) {
    vec2plane(P)
  } else {
    L <- vrotate(int, P, deg2rad(gamma))
    vec2line(L)
  }
}
