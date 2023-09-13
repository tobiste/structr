#' \deqn{R_f/\phi} 
#'
#' @param Rs numeric. 
#' @param Ri numeric. Initial or pre-deformation aspect ratio of an elliptical object
#' @param theta numeric. Fluctuation angle of undeformed marker w.r.t. computed principal direction of strain (in degree)
#'
#' @returns `phi`: the final angle (in degree) between long xis and principal strain direction, `Rf`: final axial ratio 
#' @export
Rphi <- function(Rs, Ri = 1, theta) {
  theta <- theta * pi / 180
  a <- 2 * Rs * (Ri^2 - 1) * sin(2 * theta)
  b <- (Ri^2 + 1) * (Rs - 1) + (Ri^2 - 1) * (Rs^2 + 1) * cos(2 * theta)
  tan_2phi <- a / b
  phi <- atan(tan_2phi) / 2
  
  c <- tan(phi)^2 * (1 + Ri^2 * tan(theta)^2) - Rs^2 * (tan(theta)^2 + Ri^2)
  d <- Rs^2 * tan(phi)^2 * (tan(theta)^2 + Ri^2) - (1 + Ri^2 * tan(theta)^2)
  Rf <- sqrt(c / d)
  
  return(c("phi" = phi * pi / 180, "Rf" = Rf))
}