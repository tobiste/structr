#' Simple fault analysis
#'
#' Calculates PT-axes, kinematic plane (M), and dihedra separation plane (d)
#'
#' @param x object of class `"Fault"`
#' @param ptangle numeric. angle between P and T axes in degrees (90&deg; by default).
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


  pvec <- rotate(xp, crossprod(xl, xp), -ptangle / 2)
  tvec <- rotate(xp, crossprod(xl, xp), ptangle / 2)

  p <- Line(pvec)
  t <- Line(tvec)

  mvec <- crossprod(xl, xp)
  m <- Plane(mvec)

  dvec <- crossprod(crossprod(xl, xp), xp)
  d <- Plane(dvec)

  list(p = p, t = t, m = m, d = d)
}

#' Orthogonalization of plane and line measurement
#'
#' Both the line and the plane  are rotated in opposite directions by half the angle between
#' the slip and the plane normal vector about the vector that is pperpendicular to both.
#'
#' @param x object of class `"Pair"` or `"Fault"`
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
  stopifnot(is.Pair(x))
  fvec <- Plane(x[, 1], x[, 2]) |> Vec3()
  lvec <- Line(x[, 3], x[, 4]) |> Vec3()

  misfit <- abs(pi / 2 - angle(fvec, lvec))
  # for (i in 1:length(misfit)) {
  #   if (misfit[i] > deg2rad(20)) {
  #     warning(paste("Mifit angle is", misfit[i], "degrees."))
  #   }
  # }

  # Warn if misfit > 20
  exceed_idx <- which(misfit > deg2rad(20))
  if (length(exceed_idx) > 0) {
    warning(sprintf(
      "Misfit angle exceeds 20\u00B0 for %d cases. Max misfit: %.2f\u00B0.",
      length(exceed_idx),
      rad2deg(max(misfit))
    ))
  }

  ax <- crossprod.spherical(fvec, lvec)
  ang <- (angle(lvec, fvec) - pi / 2) / 2
  list(
    fvec = rotate(fvec, ax, ang),
    lvec = rotate(lvec, ax, -ang),
    misfit = misfit
  )
}

#' @rdname pair_correct
#' @export
correct_pair <- function(x) {
  corr <- misfit_pair(x)
  p <- Plane(corr$fvec)
  l <- Line(corr$lvec)

  if (inherits(x, "fault")) {
    Fault(p[, 1], p[, 2], l[, 1], l[, 2], x[, 5], correction = FALSE)
  } else {
    Pair(p[, 1], p[, 2], l[, 1], l[, 2], correction = FALSE)
  }
}

#' Extract components of fault object
#'
#' [Fault_plane()] extracts the orientation of the fault plane,
#' [Fault_slip()] extracts the orientation of the slip vector, and
#' [Fault_rake()] extracts the rake of the fault, i.e. the angle between fault
#' slip vector and fault strike.
#'
#' @inheritParams fault_analysis
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
  stopifnot(is.Fault(x))
  strike <- Line(x[, 1] - 90, rep(0, nrow(x)))
  slip <- Line(x[, 3], x[, 4])
  x[, 5] * angle(strike, slip)
}

#' @rdname Fault_components
#' @export
Fault_slip <- function(x) {
  stopifnot(is.Pair(x))
  azi <- x[, 3]
  inc <- x[, 4]
  # sense <- x[, 5]
  # azi_cor <- ifelse(sense>1, azi + 180, azi)
  Line(azi, inc)
}

#' @rdname Fault_components
#' @export
Fault_plane <- function(x) {
  stopifnot(is.Pair(x))
  Plane(x[, 1], x[, 2])
}


#' Fault from plane and rake
#'
#' @param p object of class `"Plane"`
#' @param rake Angle in degrees in the range  −180&deg; and 180&deg;
#' @param ... optional arguments passed to [Fault()]
#'
#' @details Rake is used to describe the direction of fault motion with respect
#' to the strike (measured anticlockwise from the horizontal, up is positive; values between −180° and 180°):
#' - left-lateral strike slip: rake near 0&deg;
#' - right-lateral strike slip: rake near 180&deg;
#' - normal: rake near −90&deg;
#' - reverse/thrust: rake near +90&deg;
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' fr <- Fault_from_rake(Plane(c(120, 120, 100), c(60, 60, 50)), c(84.7202, -10, 30))
#' plot(fr, col = 1:3)
#'
#' fr2 <- Fault_from_rake(Plane(c(90, 90, 90), c(80, 40, 10)), c(10, 20, 90))
#' plot(fr2, col = 1:3)
Fault_from_rake <- function(p, rake, ...) {
  stopifnot(is.Plane(p))
  strike <- Line(p[, 1] - 90, rep(0, nrow(p)))
  rake <- rake %% 360
  rake <- ifelse(rake > 180, rake - 360, rake)

  rake_l <- rotate(strike, Line(0, 90), rake)
  l <- rotate(rake_l, strike, p[, 2])

  Fault(p[, 1], p[, 2], l[, 1], l[, 2], sense = sign(rake), ...)
}

rake2pitch <- function(rake) {
  rake <- abs(rake)
  ifelse(rake > 90, 180 - rake, rake)
}

#' Pitch of a line of a plane
#'
#' @inheritParams Fault_from_rake
#' @param l object of class `"Line"`
#'
#' @returns Pitch angle in degrees
#' @export
#'
#' @examples
#' Plane_pitch(Plane(90, 70), Line(175, 13)) # 13.9
Plane_pitch <- function(p, l) {
  Fault_rake(Fault(p[, 1], p[, 2], l[, 1], l[, 2], 1)) |> rake2pitch()
}

#' Apparent dip direction
#'
#' @inheritParams Fault_from_rake
#' @param apparent_dip angle in degrees
#'
#' @returns Azimuth in degrees
#' @export
#'
#' @examples
#' apparent_dip_direction(Plane(63, 45), 10)
apparent_dip_direction <- function(p, apparent_dip) {
  res <- Fault_from_rake(p, rake = -apparent_dip)
  res[, 3]
}

#' Apparent dip direction
#'
#' @param a1,a2 `"Line"` objects containing the apparent dips and dip directions
#'
#' @returns `"Plane"` object
#' @export
#'
#' @examples
#' a1 <- Line(45, 22)
#' a2 <- Line(352, 10)
#' res <- Plane_from_apparent_dips(a1, a2)
#'
#' stereoplot()
#' points(rbind(a1, a2))
#' lines(res, lty = 2)
Plane_from_apparent_dips <- function(a1, a2) {
  res <- crossprod(a1, a2)
  Plane(res)
}
