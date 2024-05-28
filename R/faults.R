#' Simple fault analysis
#' 
#' Calculates PT-axes, kinematic plane (M), and dihedra separation plane (d)
#'
#' @param x fault object
#' @param ptangle angle between P and T axes in degrees (90 by default) .
#'
#' @returns list
#' @export
#'
#' @examples
#' f <- Fault(120, 60, 110, 58, 1)
#' fault_analysis(f)
fault_analysis <- function(x, ptangle = 90){
 
  ptangle <- deg2rad(ptangle)
  x_corr <- correct_pair(x)  
  xp = x_corr$fvec
  xl = x_corr$lvec
  if(x[, 5] < 0) xl = -xl
  
  pvec = vrotate(xp, vcross(xl, xp), -ptangle /2)
  tvec = vrotate(xp, vcross(xl, xp), ptangle /2)
  
  p <- vec2line(pvec)
  t <- vec2line(tvec)
  
  mvec <- vcross(xl, xp)
  m <- vec2plane(mvec)
  
  dvec <- vcross(vcross(xl, xp), xp)
  d <- vec2plane(dvec)
  
  list(p = p, t = t, m = m, d = d)
}

#' Orthogonalization of plane and line measurement 
#'
#' @param x Pair or Fault
#'
#' @returns list with the orthogonalized plane and line measurements 
#' (as 3d vectors) and the misfit angles (in radians)
#' @export
#'
#' @examples
#' p <- Pair(120, 60, 110, 58)
#' correct_pair(p)
correct_pair <- function(x){
  fvec <- Plane(x[, 1], x[, 2]) |> plane2vec()
  lvec <- Line(x[, 3], x[, 4]) |> line2vec()
  
  misfit = abs(pi/2 - vangle(fvec, lvec))
  if(misfit > deg2rad(20)){
    warning(paste("Mifit angle is", misfit, "degrees."))
  }
  ax = vcross(fvec, lvec)
  ang = (vangle(lvec, fvec) - pi/2) / 2
  list(
    fvec = vrotate(fvec, ax, ang),
    lvec = vrotate(lvec, ax, -ang),
    misfit = misfit
  )
}

#' Rake of a fault
#' 
#' Angle between fault slip vector and fault strike.
#' 
#' @param x fault object
#'
#' @returns numeric. Angle in degrees
#' @export
#'
#' @examples
#' f <- Fault(120, 60, 110, 58, 1)
#' fault_rake(f)
fault_rake <- function(x){
  #x <- correct_pair(x)
  strike <- Line(x[, 1]-90, 0) |> line2vec()
  slip <- Line(x[, 3], x[, 4]) |> line2vec()
  
  rake <- x[5] * vangle(strike, slip) / DEG2RAD()
  setNames(rake, 'rake')
}
