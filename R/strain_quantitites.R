# Mohr circle for strain

#' Mohr Circle Diagram for Strain
#'
#' Plots the Mohr Circle for Strain
#'
#' @param lambda1,lambda2,lambda3 numeric. Magnitude of quadratic elongation \eqn{\lambda}
#' @param phi numeric. (optional) Angle (in degrees) for a specific strain
#' @param col color for the stress state for a given `phi`.
#' @param fg,fg12,fg23 border color for the Mohr Circles spanning lambda1-lambda3, lambda1-lambda2, and lambda2-lambda3, respectively
#' @param bg,bg12,bg23 fill color for the Mohr Circles spanning lambda1-lambda3, lambda1-lambda2, and lambda2-lambda3, respectively
#' @param full.circle logical. Should the complete Mohr circle be shown, or only
#' the upper (positive shear stress) part of the circle?
#' @param include.zero logical. the plot range be extended to include `lambda = 0`?
#' @param xlim,ylim range of plot
#' @param round integer indicating the number of decimal places to be used for rounding.  
#' @param axes logical. Show axes of plot? 
#' @param ... optional graphical parameters.
#' 
#' @returns matrix with the lambda and gamma values for given `phi`
#' 
#' @seealso [Mohr_plot()] for Stress. [strain] for converting strain quantities
#' @export
#' @examples
#' Mohr_strain(lambda1 = 4, lambda3 = 0.25, phi = 25, col = 'red')
#' (Mohr_strain(lambda1 = 4, lambda2 = 1, lambda3 = 0.25, phi = c(0, 25, 50, 45), col = 'red', full.circle = TRUE, axes = FALSE))
Mohr_strain <- function(lambda1, lambda2 = NA, lambda3,
                        #lambda_x = NA, lambda_z = NA, gamma_xz = NA,
                        phi = NULL,
                        fg = par("col"), bg = 'lightgray',
                        fg23 = par("col"), bg23 = 'white',
                        fg12 = par("col"), bg12 = 'white',
                        axes = TRUE,
                        col = "black", full.circle = FALSE, include.zero = TRUE, xlim = NULL, ylim = NULL, 
                      round = 1,
                      ...) {
  lambda_x = NA; lambda_z = NA; gamma_xz = NA
  
  phis <- seq(0, 180, length.out = 512)
  
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
      sigma1 = lambda1, sigma3 = lambda2
    )
    lambda12 <- as.numeric(stress_vec12[1, ])
    gamma12 <- as.numeric(stress_vec12[2, ])
    
    stress_vec23 <- sapply(
      X = phis, FUN = stress_transformation, 
      sigma1 = lambda2, sigma3 = lambda3
    )
    lambda23 <- as.numeric(stress_vec23[1, ])
    gamma23 <- as.numeric(stress_vec23[2, ])
    
  }
  
  ##  Expression for axes
    yLab <- bquote("Shear strain," ~ gamma*"'")
    xLab <- bquote("Quadratic elongation," ~ lambda*"'")
  
    mean_lambda <- mean(range(lambda), na.rm = TRUE)
    sigma_dev <- abs(lambda1-mean_lambda)
    
  xlim <- if (isTRUE(include.zero)) c(min(0, lambda3), lambda1) * 1.05 else xlim
  if (isFALSE(full.circle)) {
     range_f <- c(0, 1)
     } else {
       range_f <- c(-1, 1)
     }
  if(is.null(ylim)) ylim <- sigma_dev * range_f
  
  plot(
    xlim, ylim,
    type = "n",
    xlab = xLab, ylab = yLab,
    xaxs = "i",
    xlim = xlim,
    ylim = ylim,
    asp = 1,
    axes = axes
  )
  
  
  # graphics::lines(lambda, gamma, col = col, ...)
  graphics::symbols(mean_lambda, y = 0, circles = sigma_dev, inches = FALSE, fg = fg, bg = bg, add = TRUE)
  
  all_principals <- !is.na(lambda2)
  if (all_principals) {
    mean_l23 <- mean(range(lambda23))
    graphics::symbols(mean_l23, y = 0, circles = lambda2 - mean_l23, inches = FALSE, fg = fg23, bg = bg23, add = TRUE)
    # graphics::lines(lambda23, gamma23, col = col, ...)
    # graphics::lines(lambda12, gamma12, col = col, ...)
    mean_l12 <- mean(range(lambda12))
    graphics::symbols(mean_l12, y = 0, circles = lambda1 - mean_l12, inches = FALSE, fg = fg12, bg = bg12, add = TRUE)
  }
  
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
                       lty = 2, col = col)
    graphics::points(lambda_i, gamma_i, col = col)
    
    title(sub = bquote(gamma*"'"==.(round(gamma_i, round))~"|"~lambda*"'"==.(round(lambda_i, round))))
    
    return(invisible(cbind(lambda_i, gamma_i)))
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
