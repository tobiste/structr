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
#' slip vector and fault strike (measured clockwise from strike!)
#' [Fault_sense()] extracts the fault sense from the rake (1: normal, -1: reverse)
#'
#' @inheritParams fault_analysis
#' @param steps Integer. Either 2, 4, or 8 steps for parsing the fault sense
#'
#' @returns numeric. `"Plane"`, `"Ray"` or angle in degrees, respectively
#' @name Fault_components
#'
#' @examples
#' f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' Fault_plane(f)
#' Fault_slip(f)
#' Fault_rake(f)
#' Fault_sense(f, 2)
#' Fault_sense(f, 4)
#' Fault_sense(f, 8)
NULL

#' @rdname Fault_components
#' @export
Fault_rake <- function(x) {
  stopifnot(is.Fault(x))
  strike <- Line(dd2rhr(x[, 1]), rep(0, nrow(x)))
  slip <- Line(x)
  x[, 5] * angle(strike, slip)
}

#' @rdname Fault_components
#' @export
Pair_plane <- function(x) {
  stopifnot(is.Pair(x))
  # azi <- x[, 3]
  # inc <- x[, 4]
  # sense <- x[, 5]
  # # azi_cor <- ifelse(sense>1, azi + 180, azi)
  # Ray(azi, inc, sense = sense)
  Line(x)
}

#' @rdname Fault_components
#' @export
Fault_slip <- Pair_plane <- function(x) {
  stopifnot(is.Fault(x))
  # azi <- x[, 3]
  # inc <- x[, 4]
  # sense <- x[, 5]
  # # azi_cor <- ifelse(sense>1, azi + 180, azi)
  # Ray(azi, inc, sense = sense)
  Line(x)
}


#' @rdname Fault_components
#' @export
Fault_plane <- Pair_plane <- function(x) {
  stopifnot(is.Pair(x))
  Plane(x)
}

#' @rdname Fault_components
#' @export
Fault_sense <- function(x, steps = 8){
  rake <- Fault_rake(x) %% 360
  # sign(rake)
  
  if (steps == 2) {
    # Two groups: N (0–180), R (180–360)
    cut(rake,
        breaks = seq(0, 360, by = 180),
        labels = c("Normal", "Reverse"),
        include.lowest = TRUE,
        right = FALSE)
    
  } else if (steps == 4) {
    # Four groups: NS, ND, RD, RS
    cut(rake,
        breaks = seq(0, 360, by = 90),
        labels = c("Normal-sinistral", "Normal-dextral", "Reverse-dextral", "Reverse-sinistral"),
        include.lowest = TRUE,
        right = FALSE)
    
  } else if (steps == 8) {
    # Eight groups: S, NS, N, ND, D, RD, R, RS
    cut(rake,
        breaks = seq(0, 360, by = 45),
        labels = c("Sinistral", "Normal-sinistral", "Normal", "Normal-dextral", "Dextral", "Reverse-dextral", "Reverse", "Reverse-sinistral"),
        include.lowest = TRUE,
        right = FALSE)
    
  } else {
    stop("steps must be 2, 4, or 8")
  }
}



#' Fault from plane and rake
#'
#' @param p object of class `"Plane"`
#' @param rake Angle (in degrees) between fault strike and lineation. 
#' Measured clockwise from the strike, i.e. 
#' down is positive; values between 0 and 360&deg; (or −180&deg; and 180&deg;)
#' @param sense Either 1 (for normal fault movement) or -1 (reverse fault movement). 
#' Use this only if `rake` values are in range of 0 and 180&deg;. 
#' @param ... optional arguments passed to [Fault()]
#'
#' @details Rake is used to describe the direction of fault motion with respect
#' to the strike. 
#' 
#' Measured clockwise from the strike, down is positive; values between 0 and 360&deg; (or −180&deg; and 180&deg;):
#' - left-lateral strike slip: rake near 0&deg;
#' - right-lateral strike slip: rake near 180&deg;
#' - normal: rake near 90&deg;
#' - reverse/thrust: rake near -90&deg; (270&deg;)
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' fr <- Fault_from_rake(Plane(c(120, 120, 100, 0), c(60, 60, 50, 40)), 
#' rake = c(84.7202, -10, 30, 180))
#' plot(fr, col = 1:4)
#' legend("topleft", legend = Fault_sense(fr, 8), col = 1:4, pch = 16)
#'
#' fr2 <- Fault_from_rake(Plane(c(90, 90, 90), c(80, 40, 10)), 
#' rake = c(10, 20, 90), sense = c(1, 1, -1))
#' plot(fr2, col = 1:3)
#' legend("topleft", legend = Fault_sense(fr2, 8), col = 1:3, pch = 16)
Fault_from_rake <- function(p, rake, sense = NULL, ...) {
  stopifnot(is.Plane(p))
  strike <- Line(dd2rhr(p[, 1]), rep(0, nrow(p)))
  #rake_mod <- rake %% 360
  # rake <- ifelse(rake > 180, rake - 360, rake)
  if(is.null(sense)){
  qd <- (rake %% 360) <= 180
  sense = ifelse(qd, 1, -1)
  } else {
    rake <- sense * (rake %% 180)
  }
  
  l <- rotate(strike, p, rake)

  Fault(p[, 1], p[, 2], l[, 1], l[, 2], sense = sense, ...)
}

#' Fault from rake and quadrant notation
#'
#' @param p `"Plane"` object
#' @param rake numeric. Rake angle in degrees
#' @param quadrant character. Quadrant of plunge or rake direction
#' @param type character. Either `"plunge"` or `"rake"` for specifying which quadrant convention is used
#' @param sense Either 1 (for normal fault movement) or -1 (reverse fault movement). Only used when `type=="rake"`
#' 
#' @details
#' `type=="plunge"` This is the angle measured in the fault plane between the 
#' strike given by either right or left hand rule and the lineation. 
#' The angle is recorded in a clockwise sense (looking down upon the fault plane) and has a range 
#' from 0 to 180%deg;. The quadrant of plunge indicates the direction 
#' of the strike from which the angle of pitch is measured.
#' 
#' `type=="rake"` Rake is the acute angle measured in the fault plane between the strike of the fault and the 
#' lineation . Starting from the strike line, the angle is measured in a sense 
#' which is down the dip of the plane. Quadrant of pitch The variable is used to 
#' indicate the direction of the strike from which the angle of pitch is measured. 
#' Angle ranges from 0 to 90 &deg;
#' 
#' @seealso [azimuth_to_cardinal()] to convert azimuth to cardinal directions, 
#' [quadrant2dd()] to define a plane using the strike and dip quadrant notation
#'
#' @returns `"Fault"` object
#' @export
#'
#' @examples
#' dip <- c(5, 10, 15, 30, 40, 55, 65, 75, 90)
#' dip_dir <- c(180, 225, 270, 315, 360, 0, 45, 90, 135)
#' rake1 <- c(0, 45, 90, 135, 180, 45, 90, 135, 180)
#' plunge_quadrant <- c("E", "S", "W", "N", "E", "W", "E", "S", "W")
#' Fault_from_rake_quadrant(Plane(dip_dir, dip), rake1, plunge_quadrant, type = "plunge")
#' 
#' rake2 <- c(0, 45, 90, 45, 0, 45, 90, 45, 0)
#' rake_quadrant <- c("E", "S", "S", "E", "E", 'W', "N", "S", "W")
#' Fault_from_rake_quadrant(Plane(dip_dir, dip), rake2, rake_quadrant, type = "rake")
Fault_from_rake_quadrant <- function(p, rake, quadrant, type = c("plunge", "rake"), sense = NULL) {
  type <- match.arg(type)

  if (type == "plunge") {
    # try right-hand-rule
    f <- Fault_from_rake(p, rake)
    lower_hem <- f[, 4] >= 0
    to_lower <- ifelse(lower_hem, 0, 180)
    azi1 <- f[, 3] + to_lower

    # check if already correct
    res1 <- azimuth_to_cardinal(azi1, n_directions = 16)
    match1 <- sapply(seq_along(rake), function(i) grepl(quadrant[i], res1[i]))
    if (!all(match1)) {
      # do left-handrule if not
      rake2 <- rake + ifelse(match1, 0, -180)
      f <- Fault_from_rake(p, rake2)
    }
    
  } else {
    strike1 <- dd2rhr(p[, 1])
    rake_mod <- rake %% 180
    res1 <- azimuth_to_cardinal(strike1, n_directions = 16)
    match1 <- sapply(seq_along(rake_mod), function(i) grepl(quadrant[i], res1[i]))

    rake2 <- ifelse(match1, 1, -1) * rake_mod
    f <- Fault_from_rake(p, rake2 %% 180, sense = sense)
  }
  return(f)
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
