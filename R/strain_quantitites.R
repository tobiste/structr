# Mohr circle for strain

#' Mohr Circle Diagram for Strain
#'
#' Plots the Mohr Circle for Strain
#'
#' @param lambda1,lambda2,lambda3 numeric. Magnitude of quadratic elongation \eqn{\lambda}
#' @param phi numeric. (optional) Angle (in degrees) for a specific strain
#' @param col color for Mohr circle.
#' @param n integer. Resolution given amount of points along the generated path
#' representing the full Mohr circle (`512` by default).
#' @param full.circle logical. Should the complete Mohr circle be shown, or only
#' the upper (positive shear stress) part of the circle?
#' @param include.zero logical. the plot range be extended to include `lambda = 0`?
#' @param xlim,ylim range of plot
#' @param ... optional graphical parameters.
#' 
#' @seealso [Mohr_plot()] for Stress. [strain] for converting strain quantities
#' @export
#' @examples
#' Mohr_strain(lambda1 = 4, lambda3 = 0.25, phi = 25)
Mohr_strain <- function(lambda1, lambda2 = NA, lambda3,
                        #lambda_x = NA, lambda_z = NA, gamma_xz = NA,
                        phi = NULL,
                      col = "black", n = 512, full.circle = FALSE, include.zero = TRUE, xlim = NULL, ylim = NULL, 
                      digits = 1,
                      ...) {
  lambda_x = NA; lambda_z = NA; gamma_xz = NA
  
  phis <- seq(0, 180, length.out = n)
  
  stress_vec <- sapply(
    X = phis, FUN = stress_transformation, sigma_x = lambda_x, sigma_z = lambda_z,
    tau_xz = gamma_xz, sigma1 = lambda1, sigma3 = lambda3
  )
  lambda <- as.numeric(stress_vec[1, ])
  gamma <- as.numeric(stress_vec[2, ])
  
  #lambda <- (lambda1 + lambda3)/2 - (lambda3 - lambda1)/2 * cos(2*phis)
  #gamma <- (lambda3 - lambda1)/2 * sin(2*phis)
  
  if(!is.na(lambda2)){
    stress_vec12 <- sapply(
      X = phis, FUN = stress_transformation,
      sigma1 = lambda1, sigma3 = lamda2
    )
    lambda12 <- as.numeric(stress_vec12[1, ])
    gamma12 <- as.numeric(stress_vec12[2, ])
    
    stress_vec23 <- sapply(
      X = phis, FUN = stress_transformation, 
      sigma1 = lambda2, sigma3 = lamda3
    )
    lambda23 <- as.numeric(stress_vec23[1, ])
    gamma23 <- as.numeric(stress_vec23[2, ])
    
  }
  
  ##  Expression for axes
    yLab <- bquote("Shear strain," ~ gamma*"'")
    xLab <- bquote("Quadratic elongation," ~ lambda*"'")
  
  xlim <- if (include.zero) c(min(0, min(lambda, na.rm = TRUE)), max(lambda, na.rm = TRUE) * 1.05) else xlim
  ylim <- if (!full.circle) c(0, max(gamma, na.rm = TRUE)) else ylim
  
  plot(
    range(lambda), range(gamma),
    type = "n",
    xlab = xLab, ylab = yLab,
    xaxs = "i",
    xlim = xlim,
    ylim = ylim,
    asp = 1,
    axes = TRUE
  )
  
  all_principals <- !is.na(lambda2)
  
  if (all_principals) {
    graphics::lines(lambda23, gamma23, col = col, ...)
    graphics::lines(lambda12, gamma12, col = col, ...)
  }
  
 
  
  mean_lambda <- mean(range(lambda))
  
  graphics::lines(lambda, gamma, col = col, ...)
  graphics::abline(h = 0)
  
  
  if(!is.null(phi)){
    stress_veci <- sapply(
      X = 90-phi, FUN = stress_transformation, sigma_x = lambda_x, sigma_z = lambda_z,
      tau_xz = gamma_xz, sigma1 = lambda1, sigma3 = lambda3
    )
    lambda_i <- as.numeric(stress_veci[1, ])
    gamma_i <- as.numeric(stress_veci[2, ])
    
    graphics::segments(x0 = mean_lambda, y0 = 0, 
                       x1 = lambda_i, y1 = gamma_i, 
                       lty = 2)
    graphics::points(lambda_i, gamma_i)
    
    title(sub = bquote(gamma*"'"==.(round(gamma_i, digits))~"|"~lambda*"'"==.(round(lambda_i, digits))))
  }
  
  
  graphics::points(mean_lambda, 0, col = col)
}



#' Strain quantities
#' 
#' Converts different quantifications of strain.
#' 
#' @param l,l0 numeric. Final and original length, respectively
#' @param psi numeric. Angle
#' @param degree logical. Whether `psi` is given in degree (the default) or radians
#' @param e numeric. Elongation
#' @param lambda numeric. Quadratic elongation
#' 
#' @details *Longitudianal strain* (also **elongation**) is the change in length divided by the original length:
#' \deqn{e = (l-l_o)/l_o}
#' 
#' *Angular strain* is the change in angle between two lines that were initially perpendicular"
#' \deqn{\gamma = \tan \psi}
#' 
#' *Stretch* is the ratio of the final length and the original length:
#' \deqn{s = \sqrt{\lambda} = l/l_o = 1 + e}
#' 
#' *Quadratic elongation* is the quadratic stretch:
#' \deqn{\lambda = s^2 = (l/l_o)^2 = (1+e)^2}
#' 
#' @name strain
#' 
#' @examples
#' longitudinal_strain(l=10,l0=5) 
#' angular_strain(25)
#' quadratic_elongation(l=10,l0=5)
#' stretch(l=10,l0=5)
NULL

#' @rdname strain
#' @export
longitudinal_strain <- function(l, l0){
  (l - l0)/l0
}

#' @rdname strain
#' @export
angular_strain <- function(psi, degree = TRUE){
  f <- if(isTRUE(degree)) pi/180 else 1
  tan(psi * f)
}

#' @rdname strain
#' @export
quadratic_elongation <- function(e=NULL, l=NULL, l0=NULL){
  if(is.null(e)) e <- longitudinal_strain(l, l0)
  (1+e)^2
}

#' @rdname strain
#' @export
stretch <- function(e=NULL, l=NULL, l0=NULL, lambda=NULL){
  if(!is.null(e) )1 + e else if(!is.null(lambda)) sqrt(lambda)  else l/l0
}
