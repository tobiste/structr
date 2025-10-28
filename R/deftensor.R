#' Deformation Gradient Tensor
#'
#' @param Rxy,Ryz numeric. the XY and YZ strain ratio to create a strain tensor
#' with axial stretches.Values must be greater than or equal to 1.
#' @param p object of class `pair`
#' @param v1,v2 objects of class `"spherical"` or three-element vector.
#' Deformation gradient results from the rotation around axis perpendicular to
#' both vectors to rotate `v1` to `v2`.
#' @param axis,angle rotation axis and angle, axis can be an object of class
#' `"spherical"` (incl. `"line"` and `"plane"`) or a three-element vector. Angle in
#' degrees when axis is a object of class `"spherical"`, and radians otherwise.
#' @param xx,xy,xz,yx,yy,yz,zx,zy,zz numeric. Directly specify components of
#' the tensor. Identity matrix by default.
#'
#' @return 3x3 matrix.
#' @name defgrad
#'
#' @examples
#' defgrad_from_ratio(2, 3)
#' defgrad_from_axisangle(Line(120, 30), 45)
#' defgrad_from_vectors(Line(120, 30), Line(210, 60))
#' defgrad_from_pair(Pair(40, 20, 75, 16))
NULL

#' @rdname defgrad
#' @export
defgrad_from_ratio <- function(Rxy = 1, Ryz = 1) {
  stopifnot(Rxy >= 1, Ryz >= 1)
  A <- diag(3)
  y <- (Ryz / Rxy)**(1 / 3)
  A * c(y * Rxy, y, y / Ryz)
}

#' @rdname defgrad
#' @export
defgrad_from_pair <- function(p) {
  stopifnot(is.Pair(p))
  pl <- Line(p[, 3], p[, 4]) |> Vec3()
  pp <- Plane(p[, 1], p[, 2]) |> Vec3()

  D <- rbind(
    pl,
    crossprod.Vec3(pp, pl),
    pp
  )
  rownames(D) <- colnames(D) <- NULL
  t(D)
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
  )
}

#' @rdname defgrad
#' @export
defgrad_from_comp <- function(xx = 1, xy = 0, xz = 0, yx = 0, yy = 1, yz = 0,
                              zx = 0, zy = 0, zz = 1) {
  rbind(
    c(xx, xy, xz),
    c(yx, yy, yz),
    c(zx, zy, zz)
  )
}



#' Velocity gradient and Deformation gradient tensors
#'
#' Calculates the velocity gradient tensor as the matrix logarithm  of the
#' deformation gradient tensor divided by given time, and
#' the deformation gradient tensor accumulated after some time.
#'
#' @param R 3x3 matrix. Deformation gradient tensor.
#' @param V 3x3 matrix. Velocity gradient tensor.
#' @param time numeric. Total time (default is 1)
#' @param steps numeric. Time increments (default is 1)
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
#' L <- velgrad_from_defgrad(D, time = 10)
#' L
#' defgrad_from_velgrad(L, time = 10, steps = 2)
NULL

#' @rdname gradient
#' @export
velgrad_from_defgrad <- function(R, time = 1) {
  # L = pracma::logm(R) / time
  expm::logm(R) / time
}

#' @rdname gradient
#' @export
defgrad_from_velgrad <- function(V, time = 1, steps = 1) {
  if (steps > 1) {
    R <- list()
    t <- seq(0, time, steps)
    for (i in t) {
      Ri <- structure(expm::expm(V * i), class = "defgrad")
      Ri <- list(Ri)
      names(Ri) <- i
      R <- append(R, Ri)
    }
  } else {
    R <- structure(expm::expm(V * time), class = "defgrad")
  }
  R
}


#' Rate and spin of velocity gradient tensor
#'
#' @param x 3x3 matrix. Velocity gradient tensor.
#'
#' @return 3x3 matrix
#'
#' @name vel_rate
#'
#' @examples
#' R <- defgrad_from_comp(xx = 2, xy = 1, zz = 0.5)
#' L <- velgrad_from_defgrad(R, time = 10)
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
