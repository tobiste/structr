#' Stress components
#'
#' Calculates some stress components
#'
#' @param sigma1,sigma3 numeric. Magnitudes of maximum and minimum principal stress (\eqn{\sigma_1} and \eqn{\sigma_3}), respectively.
#' @param theta numeric. Angle \eqn{\theta} between fracture and \eqn{\sigma_1}.
#' @param mu numeric. Coefficient of internal friction \eqn{\mu}.
#'
#' @returns numeric
#' @name stress-comp
#'
#' @examples
#' s1 <- 1025
#' s3 <- 250
#'
#' diff_stress(s1, s3)
#' mean_stress(s1, s3)
#' shear_stress(s1, s3, theta = 35)
#' normal_stress(s1, s3, theta = 35)
#' fracture_angle(mu = 0.6)
NULL

#' @rdname stress-comp
#' @export
diff_stress <- function(sigma1, sigma3) {
  stopifnot(sigma1 >= sigma3)
  sigma1 - sigma3
}

#' @rdname stress-comp
#' @export
mean_stress <- function(sigma1, sigma3) {
  stopifnot(sigma1 >= sigma3)
  (sigma1 + sigma3) / 2
}

#' @rdname stress-comp
#' @export
shear_stress <- function(sigma1, sigma3, theta) {
  stopifnot(sigma1 >= sigma3)
  diff_stress(sigma1, sigma3) / 2 * sind(2 * theta)
}

#' @rdname stress-comp
#' @export
normal_stress <- function(sigma1, sigma3, theta) {
  stopifnot(sigma1 >= sigma3)
  mean_stress(sigma1, sigma3) - diff_stress(sigma1, sigma3) / 2 * cosd(2 * theta)
}

#' @rdname stress-comp
#' @export
fracture_angle <- function(mu) {
  (90 + atand(mu)) / 2
}


#' Principal stresses from 2D stress components
#'
#' Determines the principal stresses and their orientations from the stress
#' components \deqn{sigma_x}, \deqn{\sigma_y}, \deqn{\tau_{xy}}.
#'
#' @inheritParams Mohr_calc
#' @param sigma_y numeric. Magnitude of normal stress acting on plane facing in Y direction
#' (\deqn{\sigma_y}).
#' @param tau_xy numeric. Magnitude of shear stress acting on planes facing X and Y
#' (\deqn{\tau_{xy}}).
#'
#' @returns angle in degrees
#'
#' @references Richard J. Lisle (1999)
#' @export
#'
#' @examples
#' PR_stress(sigma_x = 80, sigma_y = 120, tau_xy = 20)
PR_stress <- function(sigma_x, sigma_y, tau_xy) {
  mean <- (sigma_x + sigma_y) / 2
  diam <- sqrt(4 * tau_xy^2 + (sigma_x - sigma_y)^2)
  radius <- diam / 2
  s1 <- mean + radius
  s2 <- mean - radius

  theta1 <- atan(tau_xy / (s1 - sigma_y)) * 180 / pi
  return(
    data.frame(sigma1 = s1, sigma2 = s2, theta = theta1)
  )
}


#' Stress transformation
#'
#' calculates the magnitudes of the normal stress and the shear stress
#'
#' @inheritParams Mohr_calc
#' @returns A two-element list containing
#' \describe{
#' \item{\code{"normal"}}{normal stress on an inclined plane}
#' \item{\code{"shear"}}{shear stress on an inclined plane}
#' }
#' @note In addition to theta, One of the following two sets of data must be entered:
#' \enumerate{
#'                \item{`"sigma_x"`, `"sigma_z"`, `"tau_xz"`}
#'                \item{`"sigma1"`, `"sigma3"`}
#'                }
#'
#' If theta is entered in conjunction with `"sigma_x"`, `"sigma_z"`, and `"tau_xz"`,
#' it is interpreted as the angle of inclination above the horizontal.  If
#' theta is entered in conjunction with the principal stresses, then it is
#'  interpreted as the angle of inclination above the major principal plane.
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' stress_transformation(sigma_x = 80, sigma_z = 120, tau_xz = 20, theta = 78)
stress_transformation <- function(theta, sigma_x = NA, sigma_z = NA, tau_xz = NA, sigma1 = NA, sigma3 = NA) {
  ##  Convert angle to radians
  theta <- deg2rad(theta)

  ##  Calculate normal and shear stresses from sigma_x, sigma_z, tau_xz
  if (is.na(sigma_x) == FALSE && is.na(sigma_z) == FALSE && is.na(tau_xz) == FALSE) {
    sigma <- (sigma_z + sigma_x) / 2 + (sigma_z - sigma_x) / 2 * cos(2 * theta) + tau_xz * sin(2 * theta)
    tau <- (sigma_z - sigma_x) / 2 * sin(2 * theta) - tau_xz * cos(2 * theta)

    ##  Calculate normal and shear stresses from principal stresses
  } else {
    if (is.na(sigma1) == FALSE && is.na(sigma3) == FALSE) {
      sigma <- (sigma1 + sigma3) / 2 + (sigma1 - sigma3) / 2 * cos(2 * theta)
      tau <- (sigma1 - sigma3) / 2 * sin(2 * theta)
    }
  }

  return(list(normal = sigma, shear = tau))
}


#' Mohr circle parameters
#'
#' calculates parameters of the Mohr circle
#
#' @param sigma_x numeric. Magnitude of normal stress acting in the horizontal direction
#' @param sigma_z numeric. Magnitude of normal stress acting in the vertical direction
#' @param tau_xz numeric. Magnitude of shear stress acting on the same plane as `"sigma_x"`
#' @param sigma1 numeric. Magnitude of major principal stress
#' @param sigma3 numeric. Magnitude of minor principal stress
#' @param theta numeric. Angles (degrees); defaults to 0-180 in increments of 1
#' @note One of the following two sets of data must be entered
#' \enumerate{
#' \item{`"sigma_x"`, `"sigma_z"`, `"tau_xz"`}
#' \item{`"sigma1"`, `"sigma3"`}
#' }
#' @returns A five-element list containing
#' \describe{
#' \item{`"mean"`}{center of Mohr circle, i.e. the mean stress}
#' \item{`"deviatoric"`}{radius of Mohr circle, i.e. the deviatoric stress}
#' \item{`"sigma"`}{vector of normal stresses for Mohr circle}
#' \item{`"tau"`}{vector of shear stresses for Mohr circle}
#' \item{`"theta"`}{vector of angles (degrees)}
#' }
#' @details
#' If theta is entered in conjunction with `"sigma_x"`, `"sigma_z"`, and `"tau_xz"`, it is interpreted
#'          as the angle of inclination above the horizontal.  If `"theta"` is entered in conjunction
#'          with the principal stresses, then it is interpreted as the angle of inclination above the
#'          major principal plane.
#' @seealso [Mohr_plot()]
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' Mohr_calc(sigma_x = 80, sigma_z = 120, tau_xz = 20)
Mohr_calc <- function(sigma_x = NA, sigma_z = NA, tau_xz = NA, sigma1 = NA, sigma3 = NA,
                      theta = seq(from = 0, to = 180, by = 1)) {
  ##  Calculate normal and shear stresses
  stress_vec <- sapply(
    X = theta, FUN = stress_transformation, sigma_x = sigma_x, sigma_z = sigma_z,
    tau_xz = tau_xz, sigma1 = sigma1, sigma3 = sigma3
  )
  sigma <- as.numeric(stress_vec[1, ])
  tau <- as.numeric(stress_vec[2, ])

  ##  Center and radius
  ##  Calculate from sigma_x, sigma_z, tau_xz
  if (is.na(sigma_x) == FALSE && is.na(sigma_z) == FALSE && is.na(tau_xz) == FALSE) {
    C <- (sigma_x + sigma_z) / 2
    R <- sqrt((((sigma_z - sigma_x) / 2)^2) + ((tau_xz)^2))

    ##  Calculate from principal stresses
  } else {
    if (is.na(sigma1) == FALSE && is.na(sigma3) == FALSE) {
      C <- (sigma1 + sigma3) / 2
      R <- (sigma1 - sigma3) / 2
    }
  }

  return(list(mean = C, deviatoric = R, normal = sigma, shear = tau, theta = theta))
}


#' Mohr Circle plot
#'
#' plots the Mohr Circle
#'
#' @inheritParams Mohr_calc
#' @param unit character. The unit used for magnitude of stress (`"MPa"` by default).
#' @param col color for Mohr circle.
#' @param n integer. Resolution given amount of points along the generated path
#' representing the full Mohr circle (`512` by default).
#' @param ... optional graphical parameters.
#' @note One of the following two sets of data must be entered
#' \enumerate{
#' \item{`"sigma_x"`, `"sigma_z"`, `"tau_xz"`}
#' \item{`"sigma1"`, `"sigma3"`}
#' }
#' @seealso [Mohr_calc()], [ggMohr()]
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' Mohr_plot(sigma_x = 80, sigma_z = 120, unit = "kPa", tau_xz = 20, col = "#B63679", lwd = 2)
#' Mohr_plot(sigma1 = 1025, sigma3 = 250, col = "#B63679", lwd = 2)
Mohr_plot <- function(sigma_x = NA, sigma_z = NA, tau_xz = NA, sigma1 = NA, sigma3 = NA,
                      unit = "MPa", col = "black", n = 512, ...) {
  ##  Calculate normal and shear stresses
  theta <- seq(from = 0, to = 180, length.out = n)
  stress_vec <- sapply(
    X = theta, FUN = stress_transformation, sigma_x = sigma_x, sigma_z = sigma_z,
    tau_xz = tau_xz, sigma1 = sigma1, sigma3 = sigma3
  )
  sigma <- as.numeric(stress_vec[1, ])
  tau <- as.numeric(stress_vec[2, ])

  ##  Expression for axes
  xLab <- bquote("Normal stress," ~ sigma[n] ~ (.(unit)))
  yLab <- bquote("Shear stress," ~ sigma[s] ~ (.(unit)))

  plot(
    range(sigma), range(tau),
    type = "n",
    xlab = xLab, ylab = yLab,
    xaxs = "i",
    xlim = c(min(0, min(sigma)), max(sigma) * 1.1),
    asp = 1
  )

  graphics::abline(h = 0)
  graphics::lines(sigma, tau, col = col, ...)
  graphics::points(mean(range(sigma)), 0, col = col)
}

#' Principle stresses
#'
#' calculates the magnitudes and directions of the principal stresses \eqn{\sigma_1} and \eqn{\sigma_3}
#
#' @inheritParams Mohr_calc
#' @returns A four-element list containing
#' \describe{
#' \item{`"sigma1"`}{magnitude of major principal stress}
#' \item{`"sigma3"`}{magnitude of minor principal stress}
#' \item{`"theta1"`}{direction of major principal stress (degrees)}
#' \item{`"theta3"`}{direction of minor principal stress (degrees)}
#'  }
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' sigma13(sigma_x = 80, sigma_z = 120, tau_xz = 20)
sigma13 <- function(sigma_x, sigma_z, tau_xz) {
  ##  Major principal stress
  sigma1 <- ((sigma_x + sigma_z) / 2) + (sqrt((((sigma_z - sigma_x) / 2)^2) + ((tau_xz)^2)))
  ##  Angle of major principal plane
  x <- (1 / (1 + ((2 * tau_xz) / (sigma_z - sigma_x))^2))
  theta1 <- acosd(sqrt(x)) / 2

  ##  Minor principal stress
  sigma3 <- ((sigma_x + sigma_z) / 2) - (sqrt((((sigma_z - sigma_x) / 2)^2) + ((tau_xz)^2)))
  ##  Angle of minor principal plane
  theta3 <- theta1 + 90

  ##  Return values
  return(list(sigma1 = sigma1, sigma3 = sigma3, theta1 = theta1, theta3 = theta3))
}



#' Maximum in-plane shear stress
#'
#' calculates the magnitude and direction of the maximum in-plane shear stress
#'
#' @inheritParams Mohr_calc
#' @returns A two-element list containing
#' \describe{
#' \item{`"tauMax"`}{maximum in-plane shear stress}
#' \item{`"theta"`}{angle of maximum in-plane shear stress (in degrees)}
#' }
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' tau_max(sigma_x = 80, sigma_z = 120, tau_xz = 20)
tau_max <- function(sigma_x, sigma_z, tau_xz) {
  ##  Maximum in-plane shear stress
  tauMax <- sqrt((((sigma_z - sigma_x) / 2)^2) + ((tau_xz)^2))
  ##  Angle of maximum in-plane shear stress
  x <- (1 / (1 + ((2 * tau_xz) / (sigma_z - sigma_x))^2))
  theta1 <- acosd(sqrt(x)) / 2
  theta <- theta1 + 45

  ##  Return values
  return(list(tau_max = tau_max, theta = theta))
}

#' Mohr Circle plot
#'
#' Plots the Mohr Circle diagram using ggplot functionality.
#'
#' @param coulomb numeric 2 element vector. Coulomb criterion containing the cohesion and the coefficient of sliding: (`c(70, 0.6)`)
#' @param sliding Sliding criteria (`0.81` by default)
#' @param units units of `sigma1`, `sigma2`, and `sigma3` (`"MPa"` by default)
#' @param sigma1,sigma2,sigma3 numeric. Magnitudes of major, intermediate, and
#' minor principal stresses. If only two principal stresses are given, only
#' one Mohr Circle will be drawn, otherwise three.
#' @param fill fill color of Mohr circle spanned by sigma1 and sigma3
#' @param alpha opacity of Mohr circle spanned sigma1 and sigma3
#' @param ... optional parameters passed to [ggforce::geom_circle()]
#'
#' @seealso [Mohr_plot()]
#'
#' @details
#' The subcaption gives the mean stress \eqn{\sigma_m}, the differential
#' stress \eqn{\sigma_d}, and the fracture angles \eqn{\theta_f} and \eqn{\alpha_f = 90 - \theta_f}
#'
#' @export
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @importFrom ggplot2 aes coord_fixed geom_abline geom_hline geom_line geom_point geom_text geom_vline ggplot labs
#' @examples
#' ggMohr(1025, 450, 250)
ggMohr <- function(sigma1, sigma2, sigma3 = NULL, coulomb = c(70, 0.6), sliding = 0.81, units = "MPa", fill = "gray", alpha = .5, ...) {
  stress <- c(sigma1, sigma2, sigma3)
  one_circle <- any(is.na(stress)) || length(stress) == 2
  stress <- sort(stats::na.omit(stress), decreasing = TRUE)

  s1 <- stress[1]
  if (one_circle) {
    s3 <- stress[2]
  } else {
    s2 <- stress[2]
    s3 <- stress[3]
  }

  circle13.r <- diff_stress(s1, s3) / 2
  circle13.m <- mean_stress(s1, s3)

  if (!one_circle) {
    circle12.r <- diff_stress(s1, s2) / 2
    circle12.m <- mean_stress(s1, s2)

    circle23.r <- diff_stress(s2, s3) / 2
    circle23.m <- mean_stress(s2, s3)
  }

  if (!is.null(coulomb)) {
    theta.f <- fracture_angle(coulomb[2]) # (90 + tectonicr:::atand(coulomb[2]))/2
  } else {
    theta.f <- 0
  }

  sigma_s <- shear_stress(s1, s3, theta.f / 2)
  sigma_n <- normal_stress(s1, s3, theta.f / 2)

  ggplot() +
    geom_circle(aes(x0 = circle13.m, y0 = 0, r = circle13.r), fill = fill, alpha = alpha, ...) +
    {
      if (!one_circle) {
        geom_circle(aes(x0 = circle23.m, y0 = 0, r = circle23.r), fill = "white", ...)
      }
    } +
    {
      if (!one_circle) {
        geom_circle(aes(x0 = circle12.m, y0 = 0, r = circle12.r), fill = "white", ...)
      }
    } +
    {
      if (!is.null(sliding)) geom_abline(intercept = 0, slope = sliding, lty = 2)
    } +
    {
      if (!is.null(coulomb)) geom_abline(intercept = coulomb[1], slope = coulomb[2], lty = 1)
    } +
    geom_point(aes(x = circle13.m, 0)) +
    geom_line(aes(x = c(circle13.m, sigma_n), y = c(0, sigma_s)), lty = 3) +
    geom_hline(yintercept = 0, alpha = .2) +
    geom_vline(xintercept = 0, alpha = .2) +
    geom_text(aes(x = (s1 + s3) / 2, y = 0), label = "sigma['m']", vjust = -.5, hjust = -1, parse = TRUE) +
    geom_text(aes(x = s3, y = 0), label = "sigma[3]", vjust = -.5, hjust = -1, parse = TRUE) +
    {
      if (!one_circle) {
        geom_text(aes(x = s2, y = 0), label = "sigma[2]", vjust = -.5, hjust = -1, parse = TRUE)
      }
    } +
    geom_text(aes(x = s1, y = 0), label = "sigma[1]", vjust = -.5, hjust = -1, parse = TRUE) +
    coord_fixed() +
    labs(
      x = bquote("Normal stress," ~ sigma[n] ~ (.(units))), y = bquote("Shear stress," ~ sigma[s] ~ (.(units))),
      caption = bquote(
        sigma[m] == .(round(circle13.m, 2)) ~ .(units) ~ "|" ~
          sigma[d] == .(round(2 * circle13.r, 2)) ~ .(units) ~ "|" ~
          theta[f] == .(round(theta.f, 2)) * degree ~ "|" ~
          alpha[f] == .(round(90 - theta.f, 2)) * degree
      )
    )
}
