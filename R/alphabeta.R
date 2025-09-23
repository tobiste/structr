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
#' @param gamma numeric. (optional) Gamma angle in degrees
#'
#' @export
#'
#' @return object of class `"Plane"`. If gamma is specified, `"Line"` object is
#' returned.
#'
#' @examples
#' # examples from Roger Marjoribanks (2016); 
#' # http://rogermarjoribanks.info/wp-content/uploads/2016/03/Plotting-alpha-to-locate-P.jpg
#' azi <- 225
#' inc <- -45
#' 
#' # single alpha-beta measurement
#' drillcore_orientation(azi, inc, 60, 320)
#' drillcore_orientation(azi, inc, 45, 220)
#'
#' # multiple alpha-beta measurements
#' my_alphas <- c(60, 45)
#' my_betas <- c(320, 220)
#' res <- drillcore_orientation(azi, inc, alpha = my_alphas, beta = my_betas)
#' 
#' 
#' # Plot core-axis, and planes in stereonet
#' plot(Line(azi, -inc), lab = "core-axis")
#' points(res, col = 2:3)
#' lines(res, col = 2:3)
#' text(res, labels = c("A", "B"), col = 2:3, pos = 4)
#' 
#' # gamma measurements
#' my_gammas <- c(0, -10)
#' res2 <- drillcore_orientation(azi, inc, alpha = my_alphas, beta = my_betas, gamma = my_gammas)
#' points(res2, col = 2:3)
drillcore_orientation <- function(azi, inc, alpha, beta, gamma = NULL) {
  # stopifnot(length(azi) == 1 & length(azi) == length(inc))
  stopifnot(length(alpha) == length(beta))
  inc <- -inc

  n <- length(alpha)
  if (length(azi) == 1) {
    azi <- rep(azi, n)
    inc <- rep(inc, n)
  }

  C <- Vec3(Line(azi, inc)) |> unclass() # core axis
  B <- Vec3(Line((azi + 180) %% 360, 90 - inc)) |> unclass() # bottom of hole line

  # E <- int <- matrix(nrow = n, ncol = 3)
  # for (i in 1:n) {
  #   if (cosd(beta[i]) <= 0) {
  #     E[i, ] <- vrotate(B[i, ], C[i, ], -deg2rad(beta[i])) |> vcross(C[i, ])
  #     int[i, ] <- vcross(E[i, ], C[i, ])
  #   } else {
  #     E[i, ] <- vrotate(B[i, ], C[i, ], deg2rad(beta[i]))
  #     int[i, ] <- vcross(C[i, ], E[i, ])
  #   }
  # }
  # P <- vrotate(C, int, -deg2rad(90 - alpha))

  # Initialize matrices
  E <- matrix(NA_real_, n, 3) # ellipse long axis
  int <- matrix(NA_real_, n, 3) # 

  beta_r <- deg2rad(beta)
  neg_cos_beta <- cos(beta_r) <= 0

  # Vectorized split computation
  # if (any(neg_cos_beta)) {
    idx_neg <- neg_cos_beta
    E[idx_neg, ] <- t(vrotate(B[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE], -beta_r[idx_neg]))
    E[idx_neg, ] <- vcross(E[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE])
    int[idx_neg, ] <- vcross(E[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE])
  # }
  # if (any(!neg_cos_beta)) {
    idx_pos <- !neg_cos_beta
    E[idx_pos, ] <- t(vrotate(B[idx_pos, , drop = FALSE], C[idx_pos, , drop = FALSE], beta_r[idx_pos]))
    int[idx_pos, ] <- vcross(C[idx_pos, , drop = FALSE], E[idx_pos, , drop = FALSE])
  # }

  P <- vrotate(C, int, -deg2rad(90 - alpha)) |>
    as.Vec3()

  if (is.null(gamma)) {
    Plane(P)
  } else {
    L <- rotate(Vec3(int), P, deg2rad(90 + gamma))
    Line(L)
  }
}
