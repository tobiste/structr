#' Simple fault analysis
#'
#' Calculates PT-axes, kinematic plane (M), and dihedra separation plane (d)
#'
#' @param x fault object
#' @param ptangle angle between P and T axes in degrees (90 by default).
#'
#' @returns list
#' @export
#'
#' @examples
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' fault_analysis(f)
fault_analysis <- function(x, ptangle = 90) {
  ptangle <- deg2rad(ptangle)
  x_corr <- misfit_pair(x)
  xp <- x_corr$fvec
  xl <- x_corr$lvec

  for (i in 1:nrow(x)) {
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
#' Both the line and the plane  are rotated in opposite directions by half the angle between 
#' the slip and the plane normal vector about the vector that is pperpendicular to both.
#'
#' @param x Pair or Fault
#'
#' @returns `misfit_pair` returns a list with the orthogonalized plane and line measurements
#' (as 3d vectors) and the misfit angles (in radians). `correct_pair` returns the orthogonalized vectors.
#' @name pair_correct
#' @examples
#' p <- Pair(120, 60, 110, 58, correction = FALSE)
#' misfit_pair(p)
#' 
#' correct_pair(p)
NULL

#' @rdname pair_correct
#' @export
misfit_pair <- function(x) {
  stopifnot(inherits(x, "pair"))
  fvec <- Plane(x[, 1], x[, 2]) |> plane2vec()
  lvec <- Line(x[, 3], x[, 4]) |> line2vec()

  misfit <- abs(pi / 2 - vangle(fvec, lvec))
  # for (i in 1:length(misfit)) {
  #   if (misfit[i] > deg2rad(20)) {
  #     warning(paste("Mifit angle is", misfit[i], "degrees."))
  #   }
  # }

  # Warn if misfit > 20
  exceed_idx <- which(misfit > deg2rad(20))
  if (length(exceed_idx) > 0) {
    warning(sprintf(
      "Misfit angle exceeds 20° for %d cases. Max misfit: %.2f°.",
      length(exceed_idx),
      rad2deg(max(misfit))
    ))
  }

  ax <- vcross(fvec, lvec)
  ang <- (vangle(lvec, fvec) - pi / 2) / 2
  list(
    fvec = vrotate(fvec, ax, ang),
    lvec = vrotate(lvec, ax, -ang),
    misfit = misfit
  )
}

#' @rdname pair_correct
#' @export
correct_pair <- function(x) {
  corr <- misfit_pair(x)
  p <- vec2plane(corr$fvec)
  l <- vec2line(corr$lvec)
  
  if(inherits(x, "fault")) {
    Fault(p[,1], p[,2], l[,1], l[,2], x[,5], correction = FALSE)
  } else {
    Pair(p[,1], p[,2], l[,1], l[,2], correction = FALSE)
  }
}

#' Extract components of fault object
#'
#' [Fault_plane()] extracts the orientation of the fault plane,
#' [Fault_slip()] extracts the orientation of the slip vector, and
#' [Fault_rake()] extracts the rake of the fault. i.e. the angle between fault slip vector and fault strike.
#'
#' @param x fault object
#'
#' @returns numeric. plane, line or angle in degrees, respectively
#' @name Fault_components
#'
#' @examples
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' Fault_plane(f)
#' Fault_slip(f)
#' Fault_rake(f)
NULL

#' @rdname Fault_components
#' @export
Fault_rake <- function(x) {
  stopifnot(inherits(x, "fault"))
  strike <- Line(x[, 1] - 90, rep(0, nrow(x)))
  slip <- Line(x[, 3], x[, 4])

  x[, 5] * vangle(strike, slip)
}

#' @rdname Fault_components
#' @export
Fault_slip <- function(x) {
  stopifnot(inherits(x, "fault"))
  azi <- x[, 3]
  inc <- x[, 4]
  # sense <- x[, 5]
  # azi_cor <- ifelse(sense>1, azi + 180, azi)
  Line(azi, inc)
}

#' @rdname Fault_components
#' @export
Fault_plane <- function(x) {
  stopifnot(inherits(x, "fault"))
  Plane(x[, 1], x[, 2])
}


#' Fault from plane and rake
#'
#' @param p plane object
#' @param rake Angle in degrees in the range  −180° and 180°
#' @param optional arguments passed to [Fault()] description
#'
#' @details Rake is used to describe the direction of fault motion with respect
#' to the strike (measured anticlockwise from the horizontal, up is positive; values between −180° and 180°):
#' - left-lateral strike slip: rake near 0°
#' - right-lateral strike slip: rake near 180°
#' - normal: rake near −90°
#' - reverse/thrust: rake near +90°
#'
#' @returns `"fault"` object
#' @export
#'
#' @examples
#' Fault_from_rake(Plane(c(120, 120, 100), c(60, 60, 50)), c(84.7202, -10, 30), ...)
Fault_from_rake <- function(p, rake) {
  strike <- Line(p[, 1] - 90, rep(0, nrow(p)))
  rake <- rake %% 360
  rake <- ifelse(rake > 180, rake - 360, rake)

  rake_l <- vrotate(strike, Line(0, 90), rake)
  l <- vrotate(rake_l, strike, p[, 2])

  Fault(p[, 1], p[, 2], l[, 1], l[, 2], sense = sign(rake), ...)
}
