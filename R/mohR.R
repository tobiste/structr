diff_stress <- function(s1, s3) {
  stopifnot(s1 >= s3)
  s1 - s3
}

mean_stress <- function(s1, s3) {
  stopifnot(s1 >= s3)
  (s1 + s3) / 2
}

shear_stress <- function(s1, s3, theta) {
  stopifnot(s1 >= s3)
  diff_stress(s1, s3) / 2 * sind(2 * theta)
}

normal_stress <- function(s1, s3, theta) {
  stopifnot(s1 >= s3)
  mean_stress(s1, s3) - diff_stress(s1, s3) / 2 * cosd(2 * theta)
}

theta <- function(phi) {
  (90 + atand(phi)) / 2
}


#' Principal stresses from 2D stress components
#' 
#' Determines the principal stresses and their orientations from the stress 
#' components \deqn{sigma_x}, \deqn{\sigma_y}, \deqn{\tau_{xy}}.
#'
#' @param sx numeric. Normal stress acting on plane facing in X direction 
#' (\deqn{sigma_x}).
#' @param sy numeric. Normal stress acting on plane facing in Y direction 
#' (\deqn{\sigma_y}).
#' @param txy numeric. Shear stress acting on planes facing X and Y 
#' (\deqn{\tau_{xy}}).
#'
#' @returns angle in degrees
#'
#' @references Richard J. Lisle (1999)
#' @export
PR_stress <- function(sx, sy, txy) {
  mean <- (sx + sy) / 2
  diam <- sqrt(4 * txy^2 + (sx - sy)^2)
  radius <- diam / 2
  s1 <- mean + radius
  s2 <- mean - radius

  theta1 <- atan(txy / (s1 - sy)) * 180 / pi
  return(
    data.frame(s1 = s1, s2 = s2, theta1 = theta1)
    )
}


#' Stress transformation
#'
#' calculates the magnitudes of the normal stress and the shear stress
#'
#' @param theta angle of inclination (degrees)
#' @param sigmaX normal stress acting in the horizontal direction
#' @param sigmaZ normal stress acting in the vertical direction
#' @param tauXZ shear stress acting on the same plane as `"sigmaX"`
#' @param sigma1 major principal stress
#' @param sigma3 minor principal stress
#' @returns A two-element list containing
#' \describe{
#' \item{\code{"sigma"}}{normal stress on an inclined plane}
#' \item{\code{"tau"}}{shear stress on an inclined plane}
#' }
#' @note In addition to theta, One of the following two sets of data must be entered:
#' \enumerate{
#'                \item{`"sigmaX"`, `"sigmaZ"`, `"tauXZ"`}
#'                \item{`"sigma1"`, `"sigma3"`}
#'                }
#'
#' If theta is entered in conjunction with `"sigmaX"`, `"sigmaZ"`, and `"tauXZ"`,
#' it is interpreted as the angle of inclination above the horizontal.  If
#' theta is entered in conjunction with the principal stresses, then it is
#'  interpreted as the angle of inclination above the major principal plane.
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' sigmaTrans(sigmaX = 80, sigmaZ = 120, tauXZ = 20, theta = 78)
sigmaTrans <- function(theta, sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA) {
  ##  Convert angle to radians
  theta <- deg2rad(theta)

  ##  Calculate normal and shear stresses from sigmaX, sigmaZ, tauXZ
  if (is.na(sigmaX) == FALSE && is.na(sigmaZ) == FALSE && is.na(tauXZ) == FALSE) {
    sigma <- (sigmaZ + sigmaX) / 2 + (sigmaZ - sigmaX) / 2 * cos(2 * theta) + tauXZ * sin(2 * theta)
    tau <- (sigmaZ - sigmaX) / 2 * sin(2 * theta) - tauXZ * cos(2 * theta)

    ##  Calculate normal and shear stresses from principal stresses
  } else {
    if (is.na(sigma1) == FALSE && is.na(sigma3) == FALSE) {
      sigma <- (sigma1 + sigma3) / 2 + (sigma1 - sigma3) / 2 * cos(2 * theta)
      tau <- (sigma1 - sigma3) / 2 * sin(2 * theta)
    }
  }

  return(list(sigma = sigma, tau = tau))
}


#' Mohr circle parameters
#'
#' calculates parameters of the Mohr circle
#
#' @param sigmaX normal stress acting in the horizontal direction
#' @param sigmaZ normal stress acting in the vertical direction
#' @param tauXZ shear stress acting on the same plane as `"sigmaX"`
#' @param sigma1 major principal stress
#' @param sigma3 minor principal stress
#' @param theta vector of angles (degrees); defaults to 0-180 in increments of 1
#' @note One of the following two sets of data must be entered
#' \enumerate{
#' \item{`"sigmaX"`, `"sigmaZ"`, `"tauXZ"`}
#' \item{`"sigma1"`, `"sigma3"`}
#' }
#' @returns A five-element list containing
#' \describe{
#' \item{`"C"`}{center of Mohr circle}
#' \item{`"R"`}{radius of Mohr circle}
#' \item{`"sigma"`}{vector of normal stresses for Mohr circle}
#' \item{`"tau"`}{vector of shear stresses for Mohr circle}
#' \item{`"theta"`}{vector of angles (degrees)}
#' }
#' @details
#' If theta is entered in conjunction with `"sigmaX"`, `"sigmaZ"`, and `"tauXZ"`, it is interpreted
#'          as the angle of inclination above the horizontal.  If `"theta"` is entered in conjunction
#'          with the principal stresses, then it is interpreted as the angle of inclination above the
#'          major principal plane.
#' @seealso [MohrCircle.plot()]
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' MohrCircle.calc(sigmaX = 80, sigmaZ = 120, tauXZ = 20)
MohrCircle.calc <- function(sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA,
                            theta = seq(from = 0, to = 180, by = 1)) {
  ##  Calculate normal and shear stresses
  stress.vec <- sapply(
    X = theta, FUN = sigmaTrans, sigmaX = sigmaX, sigmaZ = sigmaZ,
    tauXZ = tauXZ, sigma1 = sigma1, sigma3 = sigma3
  )
  sigma <- as.numeric(stress.vec[1, ])
  tau <- as.numeric(stress.vec[2, ])

  ##  Center and radius
  ##  Calculate from sigmaX, sigmaZ, tauXZ
  if (is.na(sigmaX) == FALSE && is.na(sigmaZ) == FALSE && is.na(tauXZ) == FALSE) {
    C <- (sigmaX + sigmaZ) / 2
    R <- sqrt((((sigmaZ - sigmaX) / 2)^2) + ((tauXZ)^2))

    ##  Calculate from principal stresses
  } else {
    if (is.na(sigma1) == FALSE && is.na(sigma3) == FALSE) {
      C <- (sigma1 + sigma3) / 2
      R <- (sigma1 - sigma3) / 2
    }
  }

  return(list(C = C, R = R, sigma = sigma, tau = tau, theta = theta))
}


#' Mohr Circle plot
#'
#' plots the Mohr Circle
#'
#' @param sigmaX normal stress acting in the horizontal direction
#' @param sigmaZ normal stress acting in the vertical direction
#' @param tauXZ shear stress acting on the same plane as `"sigmaX"`
#' @param sigma1 major principal stress
#' @param sigma3 minor principal stress
#' @param metric logical variable: `TRUE` for metric units (kPa), and `FALSE` for English units
#' @note One of the following two sets of data must be entered
#' \enumerate{
#' \item{`"sigmaX"`, `"sigmaZ"`, `"tauXZ"`}
#' \item{`"sigma1"`, `"sigma3"`}
#' }
#' @seealso [MohrCircle.calc()], [ggMohr()]
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' MohrCircle.plot(sigmaX = 80, sigmaZ = 120, tauXZ = 20, metric = FALSE)
MohrCircle.plot <- function(sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA,
                            metric = TRUE) {
  ##  Calculate normal and shear stresses
  theta <- seq(from = 0, to = 180, by = 1)
  stress.vec <- sapply(
    X = theta, FUN = sigmaTrans, sigmaX = sigmaX, sigmaZ = sigmaZ,
    tauXZ = tauXZ, sigma1 = sigma1, sigma3 = sigma3
  )
  sigma <- as.numeric(stress.vec[1, ])
  tau <- as.numeric(stress.vec[2, ])

  ##  Expression for axes
  if (metric == TRUE) {
    xLab <- expression("Normal Stress, " * sigma * " (kPa)")
    yLab <- expression("Shear Stress, " * tau * " (kPa)")
  } else {
    if (metric == FALSE) {
      xLab <- expression("Normal Stress, " * sigma * " (psf)")
      yLab <- expression("Shear Stress, " * tau * " (psf)")
    }
  }

  graphics::par(mgp = c(2.8, 1, 0), las = 1)
  plot(sigma, tau,
    col = "black", type = "l", xlab = xLab, ylab = yLab, xaxs = "i",
    xlim = c(min(0, min(sigma)), max(sigma) * 1.1), main = "Mohr Circle Plot", asp = 1
  )
  graphics::abline(h = 0)
}

#' Principle stresses
#'
#' calculates the magnitudes and directions of the principal stresses S1 and S2
#
#' @param sigmaX normal stress acting in the horizontal direction
#' @param sigmaZ normal stress acting in the vertical direction
#' @param tauXZ shear stress acting on the same plane as `"sigmaX"`
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
#' sigma13(sigmaX = 80, sigmaZ = 120, tauXZ = 20)
sigma13 <- function(sigmaX, sigmaZ, tauXZ) {
  ##  Major principal stress
  sigma1 <- ((sigmaX + sigmaZ) / 2) + (sqrt((((sigmaZ - sigmaX) / 2)^2) + ((tauXZ)^2)))
  ##  Angle of major principal plane
  x <- (1 / (1 + ((2 * tauXZ) / (sigmaZ - sigmaX))^2))
  theta1 <- acosd(sqrt(x)) / 2

  ##  Minor principal stress
  sigma3 <- ((sigmaX + sigmaZ) / 2) - (sqrt((((sigmaZ - sigmaX) / 2)^2) + ((tauXZ)^2)))
  ##  Angle of minor principal plane
  theta3 <- theta1 + 90

  ##  Return values
  return(list(sigma1 = sigma1, sigma3 = sigma3, theta1 = theta1, theta3 = theta3))
}



#' Maximum in-plane shear stress
#'
#' calculates the magnitude and direction of the maximum in-plane shear stress
#'
#' @param sigmaX normal stress acting in the horizontal direction
#' @param sigmaZ normal stress acting in the vertical direction
#' @param tauXZ shear stress acting on the same plane as `"sigmaX"`
#' @returns A two-element list containing
#' \describe{
#' \item{`"tauMax"`}{maximum in-plane shear stress}
#' \item{`"theta"`}{angle of maximum in-plane shear stress (in degrees)}
#' }
#' @author Kyle Elmy and Jim Kaklamanos
#' @export
#' @examples
#' tauMax(sigmaX = 80, sigmaZ = 120, tauXZ = 20)
tauMax <- function(sigmaX, sigmaZ, tauXZ) {
  ##  Maximum in-plane shear stress
  tauMax <- sqrt((((sigmaZ - sigmaX) / 2)^2) + ((tauXZ)^2))
  ##  Angle of maximum in-plane shear stress
  x <- (1 / (1 + ((2 * tauXZ) / (sigmaZ - sigmaX))^2))
  theta1 <- acosd(sqrt(x)) / 2
  theta <- theta1 + 45

  ##  Return values
  return(list(tauMax = tauMax, theta = theta))
}

#' Plot Mohr Circle in ggplots
#'
#' @param s1,s2,s3 numeric. magnitudes of sigma 1, 2, and 3, respectively.
#' @param coulomb numeric 2 element vector. Coulomb criteria (`c(70, 0.6)`)
#' @param sliding Sliding criteria (0.81 by default)
#' @param units units of the `s1`, `s2`, and `s3` (`"MPa"` by default).
#'
#' @export
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @examples
#' ggMohr(1025, 400, 250)
ggMohr <- function(s1, s2, s3, coulomb = c(70, 0.6), sliding = 0.81, units = "MPa") {
  circle13.r <- diff_stress(s1, s3) / 2
  circle13.m <- mean_stress(s1, s3)

  circle12.r <- diff_stress(s1, s2) / 2
  circle12.m <- mean_stress(s1, s2)

  circle23.r <- diff_stress(s2, s3) / 2
  circle23.m <- mean_stress(s2, s3)

  theta.f <- theta(coulomb[2]) # (90 + tectonicr:::atand(coulomb[2]))/2

  sigma_s <- shear_stress(s1, s3, theta.f / 2)
  sigma_n <- normal_stress(s1, s3, theta.f / 2)
  # theta.slope <- -atan(2*theta.f)

  ggplot2::ggplot() +
    ggforce::geom_circle(aes(x0 = circle13.m, y0 = 0, r = circle13.r), fill = "gray", alpha = .5) +
    ggforce::geom_circle(aes(x0 = circle23.m, y0 = 0, r = circle23.r), fill = "white") +
    ggforce::geom_circle(aes(x0 = circle12.m, y0 = 0, r = circle12.r), fill = "white") +
    ggplot2::geom_abline(intercept = 0, slope = sliding, lty = 2) +
    ggplot2::geom_abline(intercept = coulomb[1], slope = coulomb[2], lty = 1) +
    ggplot2::geom_point(aes(x = circle13.m, 0)) +
    ggplot2::geom_line(aes(x = c(circle13.m, sigma_n), y = c(0, sigma_s)), lty = 3) +
    ggplot2::geom_hline(yintercept = 0, alpha = .2) +
    ggplot2::geom_vline(xintercept = 0, alpha = .2) +
    ggplot2::geom_text(aes(x = s3, y = 0), label = expression(sigma[3]), vjust = -.5, hjust = -1) +
    ggplot2::geom_text(aes(x = s2, y = 0), label = expression(sigma[2]), vjust = -.5, hjust = -1) +
    ggplot2::geom_text(aes(x = s1, y = 0), label = expression(sigma[1]), vjust = -.5, hjust = -1) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = bquote(sigma[n] ~ (.(units))), y = bquote(sigma[s] ~ (.(units))), caption = bquote(theta[f] == .(round(theta.f, 2)) ~ alpha[f] == .(round(90 - theta.f, 2)))) +
    ggplot2::theme_classic()
}
