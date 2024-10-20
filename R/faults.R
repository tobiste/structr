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
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' fault_analysis(f)
fault_analysis <- function(x, ptangle = 90) {
  ptangle <- deg2rad(ptangle)
  x_corr <- correct_pair(x)
  xp <- x_corr$fvec
  xl <- x_corr$lvec
  
  for(i in 1:nrow(x)){
    if (x[i, 5] < 0) xl[i, ] <- -xl[i, ]
  }


  pvec <- vrotate(xp, vcross(xl, xp), -ptangle / 2)
  tvec <- vrotate(xp, vcross(xl, xp), ptangle / 2)

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
correct_pair <- function(x) {
  fvec <- Plane(x[, 1], x[, 2]) |> plane2vec()
  lvec <- Line(x[, 3], x[, 4]) |> line2vec()

  misfit <- abs(pi / 2 - vangle(fvec, lvec))
  for(i in 1:length(misfit)){
    if (misfit[i] > deg2rad(20)) {
      warning(paste("Mifit angle is", misfit[i], "degrees."))
    }
  }
  ax <- vcross(fvec, lvec)
  ang <- (vangle(lvec, fvec) - pi / 2) / 2
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
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' fault_rake(f)
fault_rake <- function(x) {
  # x <- correct_pair(x)
  strike <- Line(x[, 1] - 90, rep(0, nrow(x)))
  slip <- Line(x[, 3], x[, 4])

  x[, 5] * vangle(strike, slip)
  # setNames(rake, 'rake')
}


#' Fault from plane and rake
#'
#' @param p plane object
#' @param rake Angle in degrees
#' @param sense (optional) integer. Sense of the line on a fault plane. Either 1or -1 for normal or thrust offset, respectively.
#'
#' @returns returns a `"Pair"` object if `sense=NULL`, otherwise a `"Fault"` object
#' @export
#'
#' @examples
#' fault_from_rake(Plane(c(120, 120, 100), c(60, 60, 50)), c(84.7202, -10, 30), sense = c(1, -1, 1))
fault_from_rake <- function(p, rake, sense = NULL) {
  strike <- Line(p[, 1] - 90, rep(0, nrow(p)))

  rake_l <- vrotate(strike, Line(0, 90), rake)

  l <- vrotate(rake_l, strike, p[, 2])
  if (is.null(sense)) {
    Pair(p[, 1], p[, 2], l[, 1], l[, 2])
  } else {
    Fault(p[, 1], p[, 2], l[, 1], l[, 2], sense)
  }
}
