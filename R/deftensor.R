#' Deformation Gradient Tensor
#'
#' @param Rxy,Ryz numeric. the XY and YZ strain ratio to create a strain tensor
#' with axial stretches.Values must be greater than or equal to 1.
#' @param x object of class `"Pair"`, `"velgrad"` or 3x3 `"matrix"`
#' @param v1,v2 spherical objects.
#' Deformation gradient results from the rotation around axis perpendicular to
#' both vectors to rotate `v1` to `v2`.
#' @param axis,angle rotation axis and angle, axis can be an object of class `"Vec3"`,
#' `"Line"`, `"Ray"`, or `"Plane"`, or a three-element vector. Angle in
#' degrees when axis is a object of class `"Line"`, `"Ray"`, or `"Plane"`, and radians otherwise.
#' @param xx,xy,xz,yx,yy,yz,zx,zy,zz numeric. Directly specify components of
#' the tensor. Identity matrix by default.
#' @param time numeric. Total time (default is 1)
#' @param steps numeric. Time increments (default is 1)
#' @param object 3x3 `"matrix"`
#' @param ... parameters passed to function call
#'
#' @return 3x3 matrix.
#' @name defgrad
#'
#' @examples
#' defgrad_from_ratio(2, 3)
#' defgrad_from_axisangle(Line(120, 50), 60)
#' defgrad_from_vectors(Line(120, 50), Line(270, 80))
#' defgrad(Pair(40, 20, 75, 16))
NULL

#' @rdname defgrad
#' @export
is.defgrad <- function(x) inherits(x, "defgrad")

#' @rdname defgrad
#' @export
as.defgrad <- function(object) {
  stopifnot(is.matrix(object))
  structure(
    object,
    class = append(class(object), "defgrad"),
    dimnames = list(rownames(object), colnames(object))
  )
}

#' @rdname defgrad
#' @export
defgrad <- function(x, ...) UseMethod("defgrad")

#' @rdname defgrad
#' @export
defgrad_from_ratio <- function(Rxy = 1, Ryz = 1) {
  stopifnot(Rxy >= 1, Ryz >= 1)
  A <- diag(3)
  y <- (Ryz / Rxy)**(1 / 3)
  D <- A * c(y * Rxy, y, y / Ryz)
  as.defgrad(D)
}

#' @rdname defgrad
#' @export
defgrad.Pair <- function(x, ...) {
  # stopifnot(is.Pair(p))
  pl <- Line(x[, 3], x[, 4]) |> Vec3()
  pp <- Plane(x[, 1], x[, 2]) |> Vec3()

  D <- rbind(
    pl,
    crossprod.Vec3(pp, pl),
    pp
  ) |>
    unclass()
  rownames(D) <- colnames(D) <- NULL
  t(D) |> as.defgrad()
}

# defgrad_from_pairs <- function(p1, p2, symmetry = FALSE) {
#   p1l <- Line(p1[3], p1[4]) |> to_vec()
#   p1p <- Plane(p1[1], p1[2]) |> plane2vec()
#   p2l <- Line(p2[3], p2[4]) |> to_vec()
#   p2p <- Plane(p2[1], p2[2]) |> plane2vec()
#
#   R4 <- rbind(
#     Pair(p2p, p2l) %*% (p1),
#     Pair(-p2p, p2l) %*% (p1),
#     Pair(p2p, -p2l) %*% (p1),
#     Pair(-p2p, -p2l) %*% (p1)
#   )
# }

#' @rdname defgrad
#' @export
defgrad_from_vectors <- function(v1, v2) {
  v1 <- Vec3(v1)
  v2 <- Vec3(v2)

  defgrad_from_axisangle(
    axis = crossprod.Vec3(v1, v2),
    angle = angle(v1, v2) #* 180 / pi
  )
}

#' @rdname defgrad
#' @export
defgrad_from_axisangle <- function(axis, angle) {
  isvec3 <- is.Vec3(axis)
  if (!isTRUE(isvec3)) {
    axis <- Vec3(axis)
    angle <- deg2rad(angle)
  }

  x <- axis[1, 1]
  y <- axis[1, 2]
  z <- axis[1, 3]


  c <- sin(angle)
  s <- cos(angle)
  xs <- x * s
  ys <- y * s
  zs <- z * s

  xc <- x * (1 - c)
  yc <- y * (1 - c)
  zc <- z * (1 - c)
  xyc <- x * yc
  yzc <- y * zc
  zxc <- z * xc

  rbind(
    c(x * xc + c, xyc - zs, zxc + ys),
    c(xyc + zs, y * yc + c, yzc - xs),
    c(zxc - ys, yzc + xs, z * zc + c)
  ) |> 
    as.defgrad()
}

#' @rdname defgrad
#' @export
defgrad_from_comp <- function(xx = 1, xy = 0, xz = 0, yx = 0, yy = 1, yz = 0,
                              zx = 0, zy = 0, zz = 1) {
  rbind(
    c(xx, xy, xz),
    c(yx, yy, yz),
    c(zx, zy, zz)
  ) |> as.defgrad()
}

#' @rdname defgrad
#' @export
defgrad.default <- function(x, ...) {
  as.defgrad(x)
}

#' @rdname defgrad
#' @export
defgrad.velgrad <- function(x, time, steps) {
  if (steps > 1) {
    R <- list()
    t <- seq(0, time, steps)
    for (i in t) {
      Ri <- as.defgrad(expm::expm(x * i))
      Ri <- list(Ri)
      names(Ri) <- i
      R <- append(R, Ri)
    }
  } else {
    R <- as.defgrad(expm::expm(x * time))
  }
  
  return(R)
}



#' Velocity gradient and Deformation gradient tensors
#'
#' Calculates the velocity gradient tensor as the matrix logarithm  of the
#' deformation gradient tensor divided by given time, and
#' the deformation gradient tensor accumulated after some time.
#'
#' @param x 3x3 matrix. Deformation gradient tensor.
#' @param object 3x3 `"matrix"`
#' @param time numeric. Total time (default is 1)
#' @param ... parameters passed to function call
#'
#' @name gradient
#'
#' @return 3x3 matrix. If steps is > 1, then a list
#' of matrices is returned.
#'
#' @importFrom expm logm expm
#'
#' @examples
#' D <- defgrad_from_comp(xx = 2, xy = 1, zz = 0.5)
#' L <- velgrad(D, time = 10)
#' print(L)
#'
#' defgrad(L, time = 10, steps = 2)
NULL

#' @rdname gradient
#' @export
is.velgrad <- function(x) inherits(x, "velgrad")

#' @rdname gradient
#' @export
as.velgrad <- function(object) {
  stopifnot(is.matrix(object))
  structure(
    object,
    class = append(class(object), "velgrad"),
    dimnames = list(rownames(object), colnames(object))
  )
}


#' @rdname gradient
#' @export
velgrad <- function(x, ...) UseMethod('velgrad')

#' @rdname gradient
#' @export
velgrad.default <- function(x, ...) {
  as.velgrad(x)
}

#' @rdname gradient
#' @export
velgrad.defgrad <- function(x, time = 1, ...) {
  # L = pracma::logm(R) / time
  L <- expm::logm(x) / time
  as.velgrad(L)
}



#' Rate and spin of velocity gradient tensor
#'
#' @param x 3x3 matrix. Velocity gradient tensor.
#'
#' @return 3x3 matrix
#'
#' @name vel_rate
#' @seealso [velgrad()]
#'
#' @examples
#' R <- defgrad_from_comp(xx = 2, xy = 1, zz = 0.5)
#' L <- velgrad(R, time = 10)
#' velgrad_rate(L)
#' velgrad_spin(L)
NULL

#' @rdname vel_rate
#' @export
velgrad_rate <- function(x) {
  (x + t(x)) / 2
}

#' @rdname vel_rate
#' @export
velgrad_spin <- function(x) {
  (x - t(x)) / 2
}
