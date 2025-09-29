## input/output is cartesian vector (3-column vectors) and input/output angles in radians

#' Vector math operations
#'
#' @param x,y objects of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param rotaxis Axis of rotation given as object of class `"Vec3"`, `"Line"`, or `"Plane"`.
#' @param rotangle Angle of rotation in radians for `"Vec3"` objects and in degrees for `"Line"` and `"Plane"` objects.
#' @param A numeric 3x3 matrix. Transformation matrix.
#' @param norm logical. If `TRUE`, the transformed vectors are normalized to unit length.
#' @param ... arguments passed to function call
#'
#' @details
#' \describe{
#' \item{`vector_length`}{the length of a vector}
#' \item{`crossprod`}{the cross-product of two vectors, i.e. the vector perpendicular to the 2 vectors. If `y = NULL` is taken to be the sam,e vector as `x`.}
#' \item{`dotprod`}{the dot product of two vectors}
#' \item{`rotate`}{rotation of a vector about a specified vector by a specified angle}
#' \item{`angle`}{angle between two vectors}
#' \item{`project`}{projection of one vector onto the other (changes the vector
#' length of second vector, unless their are unit vectors)}
#' \item{`transform_linear`}{Linear transformation of a vector by a 3x3 matrix}
#' }
#'
#' @returns objects of same class as `x`, i.e. one of `"Vec3"`, `"Line"`, or
#' `"Plane"`. `vector_length()` and `%*%` return a real number. `angle()`
#' returns a numeric angle (in degrees, unless `x` is class `"Vec3"`).
#' @name vecmath
#'
#' @examples
#' vec1 <- Vec3(1, 0, 0)
#' vec2 <- Vec3(0, 0, 1)
#'
#' vector_length(vec1) # length of a vector
#' crossprod(vec1, vec2) # cross product
#' dotprod(vec1, vec2) # dot product
#' rotate(vec1, vec2, pi / 2) # rotation
#' angle(vec1, vec2) # angle between vectors
#' project(vec1, vec2) # projection of a vector
#' transform_linear(vec1, matrix(runif(9), 3, 3)) # linear transformation
NULL

#' @keywords internal
vlength <- function(x) {
  res <- sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2) # length of a vector
  unname(res)
}

#' @export
#' @rdname vecmath
vector_length <- function(x) {
  Vec3(x) |> vlength()
}

#' @keywords internal
vsum <- function(x) {
  cbind(x = sum(x[, 1]), y = sum(x[, 2]), z = sum(x[, 3])) # vector sum
}


#' @keywords internal
vnorm <- function(x) {
  x / vlength(x)
}

#' @keywords internal
vscale <- function(x, l) {
  l / vnorm(x) * x
}


#' @keywords internal
vcross <- function(x, y) {
  cbind(
    x = x[, 2] * y[, 3] - x[, 3] * y[, 2],
    y = x[, 3] * y[, 1] - x[, 1] * y[, 3],
    z = x[, 1] * y[, 2] - x[, 2] * y[, 1]
  )
}

# #' @rdname vecmath
# #' @export
# #' @keywords internal
# crossprod <- function(x, y, ...) UseMethod("crossprod")

# #' @export
# crossprod.default <- function(x, y, ...) base::crossprod(x, y, ...)

#' @rdname vecmath
#' @exportS3Method base::crossprod
crossprod.spherical <- function(x, y = NULL, ...) {
  xv <- Vec3(x) |> unclass()
  if (is.null(y)) yv <- xv else yv <- Vec3(y) |> unclass()

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

# `%x%.spherical` <- function(x, y) crossprod.spherical(x, y)

vdot <- function(x, y) {
  # equivalent to: x %*% t(y)
  res <- x[, 1] * y[, 1] + x[, 2] * y[, 2] + x[, 3] * y[, 3]
  unname(res)
}

# #' @exportS3Method base::`%*%`
# #' @rdname vecmath
# `%*%.spherical` <- function(x, y) {
#   dotprod(x, y)
# }

#' @export
#' @rdname vecmath
dotprod <- function(x, y) {
  e1 <- Vec3(x) |> unclass()
  e2 <- Vec3(y) |> unclass()

  vdot(e1, e2)
}

# `%*%` <- function(x, y) UseMethod("`%*%`")
# `%*%.default` <- function(x, y) x %*% y

#' @keywords internal
vrotate <- function(x, rotaxis, rotangle) {
  rotaxis <- vnorm(rotaxis)
  vax <- vcross(rotaxis, x)
  x + vax * sin(rotangle) + vcross(rotaxis, vax) * 2 * (sin(rotangle / 2))^2 # Helmut
}

#' @export
#' @rdname vecmath
rotate <- function(x, rotaxis, rotangle) UseMethod("rotate")

#' @export
rotate.default <- function(x, rotaxis, rotangle) vrotate(x, rotaxis, rotangle)

#' @export
#' @rdname vecmath
#' @method rotate spherical
rotate.spherical <- function(x, rotaxis, rotangle) {
  x3 <- Vec3(x) |> unclass()
  rotaxis3 <- Vec3(rotaxis) |> unclass()

  if (is.Line(x) | is.Plane(x)) {
    rotangle <- deg2rad(rotangle)
  }

  x_rot <- vrotate(x3, rotaxis3, rotangle) |>
    Vec3()

  if (is.Line(x)) {
    x_rot <- Line(x_rot)
  } else if (is.Plane(x)) {
    x_rot <- Plane(x_rot)
  }
  x_rot
}

#' @keywords internal
vangle <- function(x, y) {
  acos(vdot(x, y))
}

#' @export
#' @rdname vecmath
angle <- function(x, y) UseMethod("angle")

#' @export
angle.default <- function(x, y) vangle(x, y)

#' @export
angle.matrix <- function(x, y) vangle(x, y)


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

#' @keywords internal
vproject_length <- function(x, y) {
  xn <- vnorm(x)
  yn <- vnorm(y)
  vdot(xn, yn)
}

#' @export
#' @rdname vecmath
project <- function(x, y) UseMethod("project")

#' @export
project.default <- function(x, y) vproject(x, y)


#' @keywords internal
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



#' @keywords internal
vreject <- function(x, y) {
  x - vproject(x, y)
}

#' @export
#' @rdname vecmath
reject <- function(x, y) UseMethod("reject")

#' @export
#' @rdname vecmath
reject.default <- function(x, y) vreject(x, y)


#' @export
#' @method reject spherical
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

#' @keywords internal
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



#' @keywords internal
vtransform <- function(x, A, norm = FALSE) {
  stopifnot(is.matrix(A) & dim(A) == c(3, 3))

  xt <- t(A %*% t(x))
  # xt <- t(vdot(A, x))
  if (norm) xt <- vnorm(xt)
  return(xt)
}

#' @export
#' @rdname vecmath
transform_linear <- function(x, A, norm = FALSE) {
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

#' @keywords internal
v_antipode <- function(x) {
  -x
}


#' Spherical Linear Interpolation (Slerp)
#'
#' Returns the spherical linear interpolation of points between two vectors
#'
#' @param x0,x1 objects of class `"Vec3"`, `"Line"`, or `"Plane"` of the first
#' and the last points of the to be interpolated arc.
#' @param t numeric. Interpolation factor(s) (`t = [0, 1]`).
#'
#' @note For non-unit vectors the interpolation is not uniform.
#'
#' @details
#' A Slerp path is the spherical geometry equivalent of a path along a line
#' segment in the plane; a great circle is a spherical geodesic.
#'
#' @export
#'
#' @examples
#' x0 <- Line(120, 7)
#' x1 <- Line(10, 13)
#' t <- seq(0, 1, .05)
#' xslerp <- slerp(x0, x1, t)
#'
#' plot(xslerp, col = assign_col(t))
#' points(rbind(x0, x1))
slerp <- function(x0, x1, t) {
  vslerp <- sapply(t, function(i) {
    slerp1(x0, x1, t = i) |> unclass()
  }) |>
    t() |>
    rbind() |>
    as.Vec3()

  if (is.Line(x0)) {
    Line(vslerp)
  } else if (is.Plane(x0)) {
    Plane(vslerp)
  } else {
    vslerp
  }
}


slerp1 <- function(x0, x1, t) {
  stopifnot(is.Vec3(x0) | is.Line(x0) | is.Plane(x0))
  stopifnot(is.Vec3(x1) | is.Line(x1) | is.Plane(x1))
  stopifnot(
    t >= 0 | t <= 1
  )

  # if (isTRUE(na.rm)) {
  #   x0 <- x0[!rowSums(!is.finite(x0)), ]
  #   x1 <- x1[!rowSums(!is.finite(x1)), ]
  # }

  vx0 <- Vec3(x0)
  vx1 <- Vec3(x1)

  theta <- dotprod(vx0, vx1)
  sin.theta <- sin(theta)

  a <- sin(((1 - t) * theta)) / sin.theta * vx0
  b <- sin(t * theta) / sin.theta * vx1

  # return as Vec3
  return(a + b)
}
