#' Orientation of structures from drill core orientation angles
#'
#' Calculates the orientation of a plane or line from internal core angles
#' (\eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma}) of oriented drill cores
#'
#' @param azi numeric. Angle between North and the borehole trajectory projected
#' to the horizontal. The angle is measured clockwise from north and has a value
#' between 0° and 360°.
#' @param inc numeric. Acute angle between the horizontal plane and the
#' trajectory of the borehole. The value of the inclination can be between
#' −90° and 90°, where **`inc`>0° corresponds to a borehole pointing downwards**.
#' @param alpha numeric. Acute dihedral angle between the geological plane and
#' the trajectory of the borehole. The angle is restricted to be between 0° and
#' 90°, where 90° corresponds to a plane perpendicular to the borehole, i.e.
#' the trajectory of the borehole is parallel to the normal vector of the plane.
#' @param beta numeric. Angle from a reference line (line of the top of the
#' roof of the borehole profile) to the lower inflexion point of the fracture
#' trace on the borehole wall, i.e. where the perimeter of the borehole is the
#' tangent of the fracture trace. The angle is measured clockwise
#' looking in the direction of the borehole trajectory and can hence be between
#' 0° and 360°.
#' @param gamma (optional) numeric. Linear feature on a plane measured in
#' **clockwise** direction from ellipse long axis at **DOWN hole end** (positive angle).
#' If measured clockwise on a plane facing UP hole, the angle is negative.
#'
#' @name drillcore
#' @encoding UTF-8
#'
#' @references Stigsson, M., & Munier, R. (2013). Orientation uncertainty goes
#' bananas: An algorithm to visualise the uncertainty sample space on stereonets
#' for oriented objects measured in boreholes. Computers and Geosciences, 56,
#' 56–61. \doi{10.1016/j.cageo.2013.03.001}
#'
#' @return object of class `"Plane"`. If `gamma` is specified, `"Pair"` object is
#' returned.
#'
#' @examples
#' # examples from Roger Marjoribanks (2016);
#' # http://rogermarjoribanks.info/wp-content/uploads/2016/03/Plotting-alpha-to-locate-P.jpg
#' azi <- 225
#' inc <- 45
#'
#' # single alpha-beta measurement
#' drillcore_transformation(azi, inc, alpha = 60, beta = 320)
#' drillcore_transformation(azi, inc, 45, 220)
#'
#' # example from Stigsson and Munier:
#' drillcore_transformation(120, 55, 50, 270)
#'
#' # multiple alpha-beta measurements
#' my_alphas <- c(60, 45)
#' my_betas <- c(320, 220)
#' res <- drillcore_transformation(azi, inc, alpha = my_alphas, beta = my_betas)
#'
#' # Plot core-axis, and planes in stereonet
#' plot(Line(azi, inc), lab = "core-axis")
#' points(res, col = 2:3)
#' lines(res, col = 2:3)
#' text(res, labels = c("A", "B"), col = 2:3, pos = 4)
#'
#' # gamma measurements
#' my_gammas <- c(20, -10)
#' res2 <- drillcore_transformation(azi, inc, my_alphas, my_betas, my_gammas)
#' plot(res2, col = 2:3)
#' text(Line(res2), labels = c("lA", "lB"), col = 2:3, pos = 4)
NULL

# drillcore_orientation <- function(azi, inc, alpha, beta, gamma = NULL) {
#   # stopifnot(length(azi) == 1 & length(azi) == length(inc))
#   stopifnot(length(alpha) == length(beta))
#   inc <- -inc
#
#   n <- length(alpha)
#   if (length(azi) == 1) {
#     azi <- rep(azi, n)
#     inc <- rep(inc, n)
#   }
#
#   C <- Vec3(Line(azi, inc)) |> unclass() # core axis
#   B <- Vec3(Line((azi + 180) %% 360, 90 - inc)) |> unclass() # bottom of hole line
#
#   # E <- int <- matrix(nrow = n, ncol = 3)
#   # for (i in 1:n) {
#   #   if (cosd(beta[i]) <= 0) {
#   #     E[i, ] <- vrotate(B[i, ], C[i, ], -deg2rad(beta[i])) |> vcross(C[i, ])
#   #     int[i, ] <- vcross(E[i, ], C[i, ])
#   #   } else {
#   #     E[i, ] <- vrotate(B[i, ], C[i, ], deg2rad(beta[i]))
#   #     int[i, ] <- vcross(C[i, ], E[i, ])
#   #   }
#   # }
#   # P <- vrotate(C, int, -deg2rad(90 - alpha))
#
#   # Initialize matrices
#   E <- matrix(NA_real_, n, 3) # ellipse long axis
#   int <- matrix(NA_real_, n, 3) #
#
#   beta_r <- deg2rad(beta)
#   neg_cos_beta <- cos(beta_r) <= 0
#
#   # Vectorized split computation
#   # if (any(neg_cos_beta)) {
#     idx_neg <- neg_cos_beta
#     E[idx_neg, ] <- t(vrotate(B[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE], -beta_r[idx_neg]))
#     E[idx_neg, ] <- vcross(E[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE])
#     int[idx_neg, ] <- vcross(E[idx_neg, , drop = FALSE], C[idx_neg, , drop = FALSE])
#   # }
#   # if (any(!neg_cos_beta)) {
#     idx_pos <- !neg_cos_beta
#     E[idx_pos, ] <- t(vrotate(B[idx_pos, , drop = FALSE], C[idx_pos, , drop = FALSE], beta_r[idx_pos]))
#     int[idx_pos, ] <- vcross(C[idx_pos, , drop = FALSE], E[idx_pos, , drop = FALSE])
#   # }
#
#   P <- vrotate(C, int, -deg2rad(90 - alpha)) |>
#     as.Vec3()
#
#   if (is.null(gamma)) {
#     Plane(P)
#   } else {
#     L <- rotate(Vec3(int), P, deg2rad(90 + gamma))
#     Line(L)
#   }
# }

drillcore_transformation_single <- function(azi, inc, alpha, beta, gamma = NULL) {
  inc <- -inc
  stopifnot(
    any(
      length(azi) == 1,
      length(inc) == 1,
      length(alpha) == 1,
      length(beta) == 1
    )
  )


  n_BH <- c(
    cosd(beta) * cosd(alpha),
    sind(beta) * cosd(alpha),
    sind(alpha)
  )


  res <- .drillcore_helper(n_BH, azi, inc)

  Plane_normal <- Line(res["trend"] + 90, res["plunge"])
  plane <- Plane(Plane_normal)

  if (is.null(gamma)) {
    return(plane)
  } else {
    if (!is.null(gamma)) {
      stopifnot(length(gamma) == 1)

      l_BH <- n_BH |>
        t() |>
        Vec3() |>
        Line() |>
        Plane() |>
        Vec3() |>
        unclass() |>
        c()

      l_G <- .drillcore_helper(l_BH, azi, inc)
    }

    line <- Line(l_G["trend"] + 90 + 180, 90 - l_G["plunge"])

    line_rot <- rotate(line, plane, gamma)
    Pair(plane, line_rot)
  }
}

#' @rdname drillcore
#' @export
drillcore_transformation <- function(azi, inc, alpha, beta, gamma = NULL) {
  n <- length(alpha)
  if (length(azi) == 1) {
    azi <- rep(azi, n)
    inc <- rep(inc, n)
  }

  res <- sapply(seq_len(n), function(i) {
    drillcore_transformation_single(azi[i], inc[i], alpha[i], beta[i], gamma[i]) |> unclass()
  }) |> t()

  if (is.null(gamma)) as.Plane(res) else as.Pair(res)
}

.drillcore_helper <- function(n_BH, azi, inc) {
  rot_Y <- matrix(
    c(cosd(90 - inc), 0, -sind(90 - inc), 0, 1, 0, sind(90 - inc), 0, cosd(90 - inc)),
    3
  )

  rot_Z <- matrix(
    c(cosd(90 - azi), sind(90 - azi), 0, -sind(90 - azi), cosd(90 - azi), 0, 0, 0, 1),
    3
  )

  n_G <- crossprod(rot_Z, crossprod(rot_Y, n_BH)) |> t()

  temp <- n_G[1, 1] / sqrt(n_G[1, 1]^2 + n_G[1, 2]^2)

  trend <-
    if (n_G[1, 2] <= 0) {
      90 + acosd(temp)
    } else {
      90 - acosd(temp)
    }

  plunge <- asind(-n_G[1, 3])

  c(trend = trend, plunge = plunge)
}
