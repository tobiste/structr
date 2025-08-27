## input/output is cartesian vector (3-column vectors) and input/output angles in radians

#' Vector math operations
#'
#' @param x,y objects of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param rotaxis Axis of rotation given as object of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param rotangle Angle of rotation in radians for `"Vec3"` objects and in degrees for `"Line"` and `"Plane"` objects.
#' @param A numeric matrix. Transformation matrix.
#' @param norm logical. If `TRUE`, the transformed vectors are normalized to unit length.
#'
#' @returns objects of class `"Vec3"`, `"Line"`, or `"Plane"`
#' @name vecmath
#'
#' @examples
#' vec1 <- Vec3(1, 0, 0)
#' vec2 <- Vec3(0, 0, 1)
#'
#' vlength(vec1)
#' vsum(vec1)
#' vnorm(vec1)
#' crossprod(vec1, vec2)
#' vdot(vec1, vec2)
#' vec1 %*% vec2
#' rotate(vec1, vec2, pi / 2)
#' angle(vec1, vec2)
#' project(vec1, vec2)
#' reject(vec1, vec2)
NULL

vlength <- function(x) {
  res <- sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2) # length of a vector
  unname(res)
}

vsum <- function(x) {
  cbind(x = sum(x[, 1]), y = sum(x[, 2]), z = sum(x[, 3])) # vector sum
}

vnorm <- function(x) {
  x / vlength(x)
}

vscale <- function(x, l) {
  l / vnorm(x) * x
}

vcross <- function(x, y) {
  cbind(
    x = x[, 2] * y[, 3] - x[, 3] * y[, 2],
    y = x[, 3] * y[, 1] - x[, 1] * y[, 3],
    z = x[, 1] * y[, 2] - x[, 2] * y[, 1]
  )
}
crossprod <- function(x, y) UseMethod("crossprod")

#' @export
crossprod.default <- function(x, y) base::crossprod(x, y)

#' @export
#' @rdname vecmath
crossprod.spherical <- function(x, y) {
  xv <- Vec3(x) |> unclass()
  yv <- Vec3(y) |> unclass()

  xy <- vcross(xv, yv) |>
    as.Vec3()

  if (is.Line(x)) {
    Line(xy)
  } else if (is.Plane(x)) {
    Plane(xy)
  } else {
    xy
  }
}

vdot <- function(x, y) {
  # equivalent to: x %*% t(y)
  res <- x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
  unname(res)
}

#' @export
#' @rdname vecmath
`%*%.spherical` <- function(x, y) {
  e1 <- Vec3(x) |> unclass()
  e2 <- Vec3(y) |> unclass()

  vdot(e1, e2)
}



vrotate <- function(x, rotaxis, rotangle) {
  rotaxis <- vec2mat(rotaxis) |> vnorm()
  vax <- vcross(rotaxis, x)
  x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2 # Helmut
}

rotate <- function(x, rotaxis, rotangle) UseMethod("rotate")

#' @export
rotate.default <- function(x, rotaxis, rotangle) vrotate(x, rotaxis, rotangle)


#' @export
#' @rdname vecmath
rotate.spherical <- function(x, rotaxis, rotangle) {
  x3 <- Vec3(x) |> unclass()
  rotaxis3 <- Vec3(rotaxis) |> unclass()

  if (is.Line(x) | is.Plane(x)) {
    rotangle <- deg2rad(rotangle)
  }

  x_rot <- vrotate(x3, rotaxis3, rotangle)

  if (is.Line(x)) {
    x_rot <- Line(x_rot)
  } else if (is.Plane(x)) {
    x_rot <- Plane(x_rot)
  } else {
    x_rot <- as.Vec3(x_rot)
  }
  x_rot
}

vangle <- function(x, y) {
  acos(vdot(x, y))
}

angle <- function(x, y) UseMethod("angle")

#' @export
angle.default <- function(x, y) vangle(x, y)

#' @export
#' @rdname vecmath
angle.spherical <- function(x, y) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  vx <- Vec3(x) |> unclass()
  vy <- Vec3(y) |> unclass()
  res <- vangle(vx, vy)
  res <- unname(res)

  if (is.Line(x) | is.Plane(x)) {
    rad2deg(res)
  } else {
    res
  }
}

vproject_length <- function(x, y) {
  xn <- vnorm(x)
  yn <- vnorm(y)
  vdot(xn, yn)
}

project <- function(x, y) UseMethod("project")

#' @export
project.default <- function(x, y) vproject(x, y)


vproject <- function(x, y) {
  vproject_length(x, y) * y
}

#' @export
#' @rdname vecmath
project.spherical <- function(x, y) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  xv <- Vec3(x) |> unclass()
  yv <- Vec3(y) |> unclass()
  xy <- vproject(xv, yv) |> Vec3()

  if (is.Line(x)) {
    Line(xy)
  } else if (is.Plane(x)) {
    Plane(xy)
  } else {
    xy
  }
}



vreject <- function(x, y) {
  x - vproject(x, y)
}

reject <- function(x, y) UseMethod("reject")

#' @export
reject.default <- function(x, y) vreject(x,y)


#' @export
#' @rdname vecmath
reject.spherical <- function(x, y) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  xv <- Vec3(x) |> unclass()
  yv <- Vec3(y) |> unclass()
  xy <- vproject(xv, yv) |> Vec3()
  rej <- xv - xy

  if (is.Line(x)) {
    Line(rej)
  } else if (is.Plane(x)) {
    Plane(rej)
  } else {
    rej
  }
}

v_orthogonalize <- function(x, y) {
  # Early exit if already orthogonal
  if (abs(vdot(x, y)) < 1e-4) {
    return(list(x, y))
  }

  # Compute orthogonalization:
  r <- vnorm(vcross(x, y)) # rotation axis
  a <- vnorm(x) # first basis
  b <- vcross(r, a) # second basis

  # Project x and y into the rotation coordinate system:
  rv1 <- c(vdot(r, x), vdot(a, x), vdot(b, x))
  rv2 <- c(vdot(r, y), vdot(a, y), vdot(b, y))

  # Compute the correction angle:
  ang <- acos(vdot(rv1, rv2))
  ang_correction <- (pi / 2 - ang) / 2

  # Rotate in opposite directions around r:
  rw1 <- vrotate(rv1, r, -ang_correction)
  rw2 <- vrotate(rv2, r, ang_correction)

  # Transform back to NED:
  w1 <- r * rw1[1] + a * rw1[2] + b * rw1[3]
  w2 <- r * rw2[1] + a * rw2[2] + b * rw2[3]

  list(w1, w2)
}



vtransform <- function(x, A, norm = FALSE) {
  stopifnot(is.matrix(A))

  xt <- t(A %*% t(x))
  # xt <- t(vdot(A, x))
  if (norm) xt <- vnorm(xt)
  return(xt)
}


#' @export
#' @rdname vecmath
transform.spherical <- function(x, A, norm = FALSE) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  x <- Vec3(x) |> unclass()
  xt <- vtransform(x, A, norm) |> Vec3()

  if (is.Line(x)) {
    Line(xt)
  } else if (is.Plane(x)) {
    Plane(xt)
  } else {
    xt
  }
}

v_antipode <- function(x) {
  -x
}
