#' Mean orientation of a set of pairs or faults
#'
#' The Frechet (geodesic \eqn{L^2}) mean and variance of a pair of foliations and lineations
#'
#' @param x object of class `"Pair"` or `"Fault"`
#' @param group character. Symmetry group of `x`. One of
#' `"orthorhombic"` (line-in-plane symmetry, e.g. foliation-lineations, cylindrical fold
#' orientations, triaxial ellipsoid orientations, and earthquake
#' focal mechanisms, and olivine),
#' `"triclinic"` (ray-in-plane symmetry, e.g. faults with slip directions),
#' `"trigonal"` (e.g. alpha-quartz),
#' `"hexagonal"` (e.g. beta-quartz), or
#' `"trivial"` (rotations).
#' If `NULL`, the group will be
#' automatically picked based on the class of `x`.
#'
#' @returns object of class `"Pair"` or `"Fault"`, respectively
#' @name mean-pair
#'
#' @examples
#' my_fault <- Fault(c("a" = 120, "b" = 120, "c" = 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' sph_mean_pair(my_fault)
#' sph_var_pair(my_fault)
NULL

#' @rdname mean-pair
#' @export
sph_mean_pair <- function(x, group = NULL) {
  res <- .sph_meanvar_pair(x, group)
  rot2pair(res$mean, fault = inherits(x, "Fault"))
}

#' @rdname mean-pair
#' @export
sph_var_pair <- function(x, group = NULL) {
  .sph_meanvar_pair(x, group)$variance
}

.sph_meanvar_pair <- function(x, group) {
  xvec <- pair2rot(x)

  if (is.null(group)) {
    group <- if (inherits(x, "Fault")) "triclinic" else "orthorhombic"
  }

  group_mat <- switch(group,
    "orthorhombic" = oriLineInPlaneGroup(),
    "triclinic" = oriRayInPlaneGroup(),
    "trivial" = oriTrivialGroup(),
    "trigonal" = oriTrigonalTrapezohedralGroup(),
    "hexagonal" = oriHexagonalTrapezohedralGroup()
  )
  ori_mean_variance(xvec, group = group_mat)
}


#' Orientation matrix from fault planes and slip directions
#'
#' Converts a set of planes and lines into a list of rotation matrices.
#'
#' @param x object of class `"Pair"` or `"Fault"`
#' @export
#'
#' @returns list of rotation matrices
#'
#' @examples
#' my_fault <- Fault(c("a" = 120, "b" = 120, "c" = 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' pair2rot(my_fault)
pair2rot <- function(x) {
  stopifnot(is.Pair(x))
  p <- Fault_plane(x) |> Vec3()
  l <- Fault_slip(x) |> Vec3()

  if (inherits(x, "Fault")) {
    l <- x[, 5] * l
  }
  cross <- crossprod(p, l)

  rotm <- rot_projected_matrix(list(pole = p, direction = l, cross = cross))
  class(rotm) <- append(class(rotm), "Rotation")
  return(rotm)
}

is.Rotation <- function(x) {
  inherits(x, "Rotation")
}

rot2pair <- function(x, fault = FALSE) {
  stopifnot(is.Rotation(x))
  mp <- as.Vec3(x[1, ])
  ml <- as.Vec3(x[2, ])

  if (fault) {
    Fault(Plane(mp), Line(ml), sense = as.integer(sign(ml[, 3])))
  } else {
    Pair(Plane(mp), Line(ml))
  }
}


rot_projected_matrix <- function(x) {
  stopifnot(is.list(x))
  n <- nrow(x$pole)

  lapply(seq_len(n), function(i) {
    m <- rbind(pole = x$pole[i, ], direction = x$direction[i, ], x$cross[i, ]) |> unclass()

    valsvecs <- eigen(crossprod(m, m), symmetric = TRUE)
    q <- valsvecs$vectors
    dSqrtInv <- valsvecs$values^(-1 / 2)
    # projection = M Q D^(-1 / 2) Q^T.
    m %*% q %*% diag(dSqrtInv) %*% t(q)
  })
}


#' The Frechet (geodesic L^2) variance of a set of orientations.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Often the mean of the rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative). The variance of the points G R about the point G center.
#' @noRd
#' @source geologyGeometry
oriVariance <- function(rs, center, group) {
  dists <- sapply(rs, oriDistance, center, group)
  sum(dists^2) / (2 * length(rs))
}


#' The Frechet (geodesic L^2) mean of a set of orientations as points in SO(3) / G.
#'
#' An iterative algorithm for computing the geodesic mean --- the orientation that minimizes the geodesic variance. The iterations continue until change-squared of epsilon is achieved or numSteps iterations have been used.
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). How many seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $mean (a special orthogonal real 3x3 matrix), $variance (a real number), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
#'
#' @noRd
#' @source geologyGeometry
#'
#' @examples
#' my_fault <- Fault(c("a" = 120, "b" = 120, "c" = 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#' my_fault_vec <- pair2rot(my_fault)
#' res1 <- ori_mean_variance(my_fault_vec, group = oriLineInPlaneGroup())
#' rot2pair(res1$mean, fault = TRUE)
#'
#' res2 <- ori_mean_variance(my_fault_vec, group = oriRayInPlaneGroup())
#' rot2pair(res2$mean, fault = TRUE)
ori_mean_variance <- function(rs, group, numSeeds = 5, numSteps = 1000, eps = sqrt(.Machine$double.eps)) {
  n <- length(rs)
  g_len <- length(group)

  # random seeds
  seeds <- if (n >= numSeeds) sample(rs, numSeeds) else sample(rs, numSeeds, replace = TRUE)

  # upper bound for variance
  best_var <- 5
  best <- vector("list", 4)

  for (seed in seeds) {
    rBar <- seed
    changeSquared <- eps + 1.0
    k <- 0

    while (changeSquared >= eps && k < numSteps) {
      w <- matrix(0, 3, 3)

      for (r in rs) {
        # compute all crossproducts in one pass
        scores <- numeric(g_len)
        mats <- vector("list", g_len)
        for (j in seq_len(g_len)) {
          cp <- crossprod(rBar, group[[j]] %*% r)
          mats[[j]] <- cp
          scores[j] <- tr(cp) # precompute trace
        }
        i <- which.max(scores)
        w <- w + rotLog(mats[[i]])
      }

      w <- w / n
      rBar <- rBar %*% rotExp(w)
      changeSquared <- tr(crossprod(w, w))
      k <- k + 1
    }

    var <- oriVariance(rs, rBar, group)
    if (var < best_var) {
      best_var <- var
      best <- list(rBar = rBar, changeSquared = changeSquared, k = k)
    }
  }

  # assign class once
  class(best$rBar) <- c("Rotation", class(best$rBar))

  list(
    variance = best_var,
    mean = best$rBar,
    changeSquared = best$changeSquared,
    numSteps = best$k
  )
}

# Previous version, kept for reference#
# oriMeanVariance <- function(rs, group, numSeeds=5, numSteps=1000, eps = sqrt(.Machine$double.eps)) {
#   seeds <- sample(rs, numSeeds, replace = (length(rs) < numSeeds))
#   # No variance is ever larger than pi^2 / 2 < 5.
#   best <- list(5)
#   for (seed in seeds) {
#     rBar <- seed
#     changeSquared <- eps + 1.0
#     k <- 0
#     while (changeSquared >= eps && k < numSteps) {
#       w <- diag(c(0, 0, 0))
#       for (r in rs) {
#         rBarTGRs <- lapply(group, function(g) crossprod(rBar, g %*% r))
#         i <- which.max(sapply(rBarTGRs, tr))
#         w <- w + rotLog(rBarTGRs[[i]])
#       }
#       w <- w / length(rs)
#       rBar <- rBar %*% rotExp(w)
#       changeSquared <- tr(crossprod(w, w))
#       k <- k + 1
#     }
#     var <- oriVariance(rs, rBar, group)
#     if (var < best[[1]])
#       best <- list(var, rBar, changeSquared, k)
#   }
#
#   class(best[[2]]) <- append(class(best[[2]]), "Rotation")
#
#   list(
#     variance=best[[1]], mean=best[[2]], changeSquared=best[[3]],
#     numSteps=best[[4]])
# }




#' The distance between two orientations as points in SO(3) / G.
#'
#' @param r,q A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @returns A real number (in the interval `[0, pi]`). The distance from G R to G R.
#' @noRd
#' @source geologyGeometry
oriDistance <- function(r, q, group) {
  qRT <- tcrossprod(q, r)
  trGQRT <- max(sapply(group, function(g) tr(g %*% qRT)))
  arcCos((trGQRT - 1) / 2)
}

#' Diameter of a set of orientations.
#'
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @returns A real number (non-negative).
#' @noRd
#' @source geologyGeometry
oriDiameter <- function(rs, group) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) oriDistance(rs[[i]], rs[[j]], group)))
  }
  max(sapply(2:(length(rs)), f))
}

#' Selecting an orientation representative near a given rotation.
#'
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A rotation matrix. The element of G R that is closest to center.
#' @noRd
#' @source geologyGeometry
oriNearestRepresentative <- function(r, center, group) {
  gRs <- lapply(group, function(g) g %*% r)
  trCTGRs <- sapply(gRs, function(gr) tr(crossprod(center, gr)))
  i <- which.max(trCTGRs)
  gRs[[i]]
}

#' Selecting orientation representatives near a given rotation.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. For each R, the element of G R that is closest to center.
#' @noRd
#' @source geologyGeometry
oriNearestRepresentatives <- function(rs, center = rs[[1]], group) {
  lapply(rs, oriNearestRepresentative, center, group)
}


#' Matrix logarithm to produce infinitesimal from finite rotation.
#'
#' @param r A 3x3 rotation matrix.
#' @return A 3x3 real matrix (anti-symmetric). The principal logarithm.
#' @noRd
#' @source geologyGeometry
rotLog <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[[4]] * rotAntisymmetricFromVector(ua[1:3])
}

#' Matrix exponentiation of infinitesimal rotation to finite rotation.
#'
#' @param w A 3x3 real matrix (anti-symmetric).
#' @return A rotation matrix.
#' @noRd
#' @source geologyGeometry
rotExp <- function(w) {
  th <- sqrt(tr(crossprod(w, w)) / 2)
  if (th == 0) {
    diag(c(1, 1, 1))
  } else {
    diag(c(1, 1, 1)) + (sin(th) / th) * w + ((1 - cos(th)) / th^2) * (w %*% w)
  }
}

# Safe sqrt, for times when you know that you don't care about arguments slightly outside [0, infinity).
squareRoot <- function(x) {
  sqrt(max(0, x))
}


#' Conversion from matrix to angle-axis representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians (between 0 and pi).
#' @noRd
#' @source geologyGeometry
rotAxisAngleFromMatrix <- function(r) {
  cosine <- (tr(r) - 1) / 2
  if (cosine >= 1) {
    # The angle is 0 and the axis doesn't matter.
    c(0, 0, -1, 0)
  } else {
    # These are 2 * sin(a) >= 0 times the true values.
    u1 <- r[3, 2] - r[2, 3]
    u2 <- r[1, 3] - r[3, 1]
    u3 <- r[2, 1] - r[1, 2]
    normU <- squareRoot(u1^2 + u2^2 + u3^2)
    if (normU != 0) {
      c(u1 / normU, u2 / normU, u3 / normU, arcCos(cosine))
    } else if (cosine > 0) {
      # The angle is 0 and the axis doesn't matter.
      c(0, 0, -1, 0)
    } else {
      # The angle is pi and the axis is unstable. So we guess the absolute
      # values of u1, u2, u3.
      u1 <- squareRoot((r[1, 1] + 1) / 2)
      u2 <- squareRoot((r[2, 2] + 1) / 2)
      u3 <- squareRoot((r[3, 3] + 1) / 2)
      # We can choose signs based on Rij = 2 ui uj when i != j.
      if (u1 != 0) {
        # Pick u1 > 0.
        if (r[1, 2] < 0) {
          u2 <- -u2
        }
        if (r[1, 3] < 0) {
          u3 <- -u3
        }
      } else if (u2 != 0) {
        # Pick u2 > 0.
        if (r[2, 1] < 0) {
          u1 <- -u1
        }
        if (r[2, 3] < 0) {
          u3 <- -u3
        }
      } else {
        # Pick u3 > 0.
        if (r[3, 1] < 0) {
          u1 <- -u1
        }
        if (r[3, 2] < 0) {
          u2 <- -u2
        }
      }
      c(u1, u2, u3, pi)
    }
  }
}

#' Conversion from vector representation to anti-symmetric matrix (infinitesimal rotation).
#'
#' @param v A real 3D vector.
#' @return A infinitsimal rotation matrix.
#' @noRd
#' @source geologyGeometry
rotAntisymmetricFromVector <- function(v) {
  matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}




### SYMMETRY GROUPS ###

# Frequently orientations are subject to some symmetry group G, which is a
# finite set of rotations (satisfying certain properties). For any rotation Q
# in G, the rotation matrix Q %*% R represents the same orientation as R does.
# Here are some common symmetry groups.

#' When trivial symmetry is used, the orientations are simply rotations. This case is so important that we have separate code, in rot.R, for doing it. But let's include the trivial group, for completeness.
oriTrivialGroup <- function() list(diag(c(1, 1, 1)))

#' Ray-in-plane symmetry is applicable to faults-with-slip-directions, such as slickensides (and certain minerals).
oriRayInPlaneGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    diag(c(-1, 1, -1))
  )
}
#' Line-in-plane symmetry is applicable to foliation-lineations, cylindrical fold orientations, triaxial ellipsoid orientations, and earthquake focal mechanisms (and olivine).
oriLineInPlaneGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    diag(c(1, -1, -1)),
    diag(c(-1, 1, -1)),
    diag(c(-1, -1, 1))
  )
}
#' Trigonal trapezohedral is the point group of alpha-quartz.
oriTrigonalTrapezohedralGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    rotMatrixAboutZ(pi * 2 / 3),
    rotMatrixAboutZ(pi * 4 / 3),
    diag(c(1, -1, -1)),
    rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 2 / 3)),
    rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 4 / 3))
  )
}
#' Hexagonal trapezohedral is the point group of beta-quartz.
oriHexagonalTrapezohedralGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    rotMatrixAboutZ(pi * 1 / 3),
    rotMatrixAboutZ(pi * 2 / 3),
    diag(c(-1, -1, 1)),
    rotMatrixAboutZ(pi * 4 / 3),
    rotMatrixAboutZ(pi * 5 / 3),
    diag(c(1, -1, -1)),
    rotMatrixAboutZ(pi * 1 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 1 / 3)),
    rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 2 / 3)),
    diag(c(-1, -1, 1)) %*% diag(c(1, -1, -1)) %*% diag(c(-1, -1, 1)),
    rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 4 / 3)),
    rotMatrixAboutZ(pi * 5 / 3) %*% diag(c(1, -1, -1)) %*%
      t(rotMatrixAboutZ(pi * 5 / 3))
  )
}
#' Replicating a list of orientations into their multiple rotation representatives.
#'
#' @param rs A list of rotation matrices. The representative rotations Rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. The set G Rs. For each orientation G R, its |G| representatives are spread out in this list.
#' @noRd
#' @source geologyGeometry
oriSymmetrizedRotations <- function(rs, group) {
  unlist(
    lapply(
      group,
      function(g) lapply(rs, function(r) g %*% r)
    ),
    recursive = FALSE, use.names = FALSE
  )
}

#' Rotation matrix about the x-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
#' @noRd
#' @source geologyGeometry
rotMatrixAboutX <- function(a) {
  cosine <- cos(a)
  sine <- sin(a)
  matrix(c(1, 0, 0, 0, cosine, sine, 0, -sine, cosine), 3, 3)
}

#' Rotation matrix about the y-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
#' @noRd
#' @source geologyGeometry
rotMatrixAboutY <- function(a) {
  cosine <- cos(a)
  sine <- sin(a)
  matrix(c(cosine, 0, -sine, 0, 1, 0, sine, 0, cosine), 3, 3)
}

#' Rotation matrix about the z-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
#' @noRd
#' @source geologyGeometry
rotMatrixAboutZ <- function(a) {
  cosine <- cos(a)
  sine <- sin(a)
  matrix(c(cosine, sine, 0, -sine, cosine, 0, 0, 0, 1), 3, 3)
}

#' Trace of a matrix.
#'
#' @param m A real nxn matrix.
#' @return A real number.
#' @noRd
#' @source geologyGeometry
tr <- function(m) {
  sum(diag(m))
}

# Safe acos, for times when you know that you don't care about arguments slightly outside [-1, 1].
#' @noRd
#' @source geologyGeometry
arcCos <- function(x) {
  acos(max(-1, min(1, x)))
}


#' Sanitizing an orientation of a fault-with-slip.
#'
#' The issue here is that geologists like describing fault slips in terms of the movement direction of the hanging wall, but the hanging wall is undefined for vertical faults. A much better-behaved convention is to describe the direction of slip using a vorticity vector. Inverse to geoPoleHangingFromPoleVorticity.
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole to a fault plane. The second row is the movement direction of the hanging wall. (The third row is the cross product of the other two.)
#' @return A real 3x3 matrix (special orthogonal). The first row is the downward-pointing pole to the fault. The second row is the vorticity vector of the slip. (The third row is the cross product of the other two.)
#' @noRd
#' @source geologyGeometry
geoPoleVorticityFromPoleHanging <- function(rot) {
  if (rot[[1, 3]] < 0) {
    rbind(rot[1, ], -rot[3, ], rot[2, ])
  } else {
    rbind(-rot[1, ], rot[3, ], rot[2, ])
  }
}

#' Desanitizing an orientation of a fault-with-slip.
#'
#' Inverse to geoPoleVorticityFromPoleHanging.
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole to a fault plane. The second row is the vorticity vector of the fault slip. (The third row is the cross product of the other two.)
#' @return A real 3x3 matrix (special orthogonal). The first row is the downward-pointing pole to the fault. The second row is the movement direction of the hanging wall. (The third row is the cross product of the other two.)
#' @noRd
#' @source geologyGeometry
geoPoleHangingFromPoleVorticity <- function(rot) {
  if (rot[[1, 3]] < 0) {
    rbind(rot[1, ], rot[3, ], -rot[2, ])
  } else {
    rbind(-rot[1, ], -rot[3, ], -rot[2, ])
  }
}

#' Trace of a matrix.
#'
#' @param m A real nxn matrix.
#' @return A real number.
#' @noRd
#' @source geologyGeometry
tr <- function(m) {
  sum(diag(m))
}

#' Cross product of three-dimensional Cartesian vectors.
#'
#' @param v A 3D real vector.
#' @param w A 3D real vector.
#' @return A 3D real vector.
#' @noRd
#' @source geologyGeometry
cross <- function(v, w) {
  c(
    v[[2]] * w[[3]] - v[[3]] * w[[2]],
    v[[3]] * w[[1]] - v[[1]] * w[[3]],
    v[[1]] * w[[2]] - v[[2]] * w[[1]]
  )
}

#' Dot product of n-dimensional Cartesian vectors.
#'
#' @param v A real vector.
#' @param w A real vector of the same dimension.
#' @return A real number.
#' @noRd
#' @source geologyGeometry
dot <- function(v, w) {
  sum(v * w)
}

#' Antipodal map sending upper-hemisphere rays to the lower hemisphere.
#'
#' In other words, negates any vector with positive `[[3]]`-component. Used in lower-hemisphere plots.
#'
#' @param v A real 3D vector.
#' @return A real 3D vector.
#' @noRd
#' @source geologyGeometry
lower <- function(v) {
  if (v[[3]] <= 0) {
    v
  } else {
    -v
  }
}
