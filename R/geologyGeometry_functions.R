# Helpers ----------------------------------------------------------------------


#' Trace of a matrix.
#'
#' @param m A real nxn matrix.
#' @return A real number.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
tr <- function(m) sum(diag(m))

# Safe acos, for times when you know that you don't care about arguments slightly outside [-1, 1].
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
arcCos <- function(x) acos(max(-1, min(1, x)))

#' Cross product of three-dimensional Cartesian vectors.
#'
#' @param v A 3D real vector.
#' @param w A 3D real vector.
#' @return A 3D real vector.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
dot <- function(v, w) sum(v * w)


#' Antipodal map sending upper-hemisphere rays to the lower hemisphere.
#'
#' In other words, negates any vector with positive `[[3]]`-component. Used in lower-hemisphere plots.
#'
#' @param v A real 3D vector.
#' @return A real 3D vector.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lower <- function(v) {
  if (v[[3]] <= 0) {
    v
  } else {
    -v
  }
}

# Safe sqrt, for times when you know that you don't care about arguments slightly outside [0, infinity).
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
squareRoot <- function(x) sqrt(max(0, x))



# Coordinate transformation ----------------------------------------------------

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
cartesianFromSpherical <- function(rpt) {
  sinPhi <- sin(rpt[[2]])
  x <- rpt[[1]] * sinPhi * cos(rpt[[3]])
  y <- rpt[[1]] * sinPhi * sin(rpt[[3]])
  z <- rpt[[1]] * cos(rpt[[2]])
  c(x, y, z)
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
sphericalFromCartesian <- function(xyz) {
  rho <- sqrt(dot(xyz, xyz))
  if (rho == 0) {
    c(0, 0, 0)
  } else {
    phi <- acos(xyz[[3]] / rho)
    rhoSinPhi <- rho * sin(phi)
    if (rhoSinPhi == 0) {
      if (xyz[[3]] >= 0) {
        c(xyz[[3]], 0, 0)
      } else {
        c(-xyz[[3]], pi, 0)
      }
    } else {
      c(rho, phi, atan2(xyz[[2]], xyz[[1]]))
    }
  }
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
cartesianFromHorizontal <- function(tz) {
  root <- sqrt(1 - tz[[2]]^2)
  c(cos(tz[[1]]) * root, sin(tz[[1]]) * root, tz[[2]])
}



# Basics -----------------------------------------------------------------------

#' Vector normalization (to have length 1).
#'
#' @param v A d-dimensional vector. Cannot be of length 0.
#'
#' @return A d-dimensional vector of length 1.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayNormalized <- function(v) v / sqrt(sum(v^2))
# similar to vnorm


#' The exponential map on the unit sphere.
#'
#' Regards `v` as a vector in the tangent space of the unit sphere at the point `p`.
#' Returns a point `q` on the unit sphere, a distance of `|v|` from `p`, in the
#' direction of v. Partial inverse to `rayLog()`.
#'
#' @param p A ray.
#' @param v A 3-dimensional vector, perpendicular to `p`, not necessarily unit.
#'
#' @return A ray.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayExp <- function(p, v) {
  normV <- sqrt(dot(v, v))
  if (normV == 0) {
    p
  } else {
    r <- rbind(p, v / normV, cross(p, v / normV))
    as.numeric(t(r) %*% c(cos(normV), sin(normV), 0))
  }
}

#' Inverse exponential map on the unit sphere.
#'
#' Returns a vector v in the tangent space to the unit sphere at p, such that
#' `rayExp(p, v) = q`.
#'
#' @param p,q Rays
#'
#' @return A 3-dimensional vector, perpendicular to p, not necessarily unit.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayLog <- function(p, q) {
  normV <- arcCos(dot(p, q))
  w <- rayNormalized(cross(p, q))
  v <- normV * cross(w, p)
  if (dot(v, q) >= 0) {
    v
  } else {
    -v
  }
}


#' Matrix logarithm to produce infinitesimal from finite rotation.
#'
#' @param r A 3x3 rotation matrix.
#' @return A 3x3 real matrix (anti-symmetric). The principal logarithm.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rotLog <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[[4]] * rotAntisymmetricFromVector(ua[1:3])
}






#' Matrix exponentiation of infinitesimal rotation to finite rotation.
#'
#' @param w A 3x3 real matrix (anti-symmetric).
#' @return A rotation matrix.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rotExp <- function(w) {
  th <- sqrt(tr(crossprod(w, w)) / 2)
  if (th == 0) {
    diag(c(1, 1, 1))
  } else {
    diag(c(1, 1, 1)) + (sin(th) / th) * w + ((1 - cos(th)) / th^2) * (w %*% w)
  }
}


#' Projection of a vector onto the plane perpendicular to another vector.
#'
#' For example, when v = c(0, 0, 1), returns the upward-most ray perpendicular
#' to pole.
#'
#' @param pole A d-dimensional vector perpendicular to the plane of interest.
#' Must have non-zero length. Need not be unit.
#'
#' @param v A d-dimensional vector. Should not be parallel to pole. Need not be
#' unit.
#'
#' @return A unit d-dimensional vector.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayOrthogonalProjection <- function(pole, v) rayNormalized(v - pole * dot(v, pole) / dot(pole, pole))

#' A ray perpendicular to a given vector, deterministically.
#'
#' The result is deterministic but arbitrary. Due to the hairy ball theorem, the
#' result cannot be chosen to depend smoothly on the input.
#'
#' @param v A 3-dimensional vector. Need not be unit.
#'
#' @return A ray, perpendicular to v.
#'
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayOrthogonal <- function(v) {
  if (abs(v[1]) < 0.5) {
    w <- c(1, 0, 0)
  } else if (abs(v[2]) < 0.5) {
    w <- c(0, 1, 0)
  } else {
    w <- c(0, 0, 1)
  }
  rayOrthogonalProjection(v, w)
}

#' A ray perpendicular to a given vector, probabilistically.
#'
#' In theory, the returned ray is uniformly chosen on the circle's worth of
#' rays perpendicular to the given vector.
#'
#' @param v A 3-dimensional vector. Need not be unit.
#'
#' @return A ray, perpendicular to v.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayOrthogonalUniform <- function(v) rayOrthogonalProjection(v, rayUniform())

#' Uniformly random points on the unit sphere.
#'
#' @param n A real number (positive integer) or NULL.
#'
#' @return If n is NULL, then a single ray. If n is a positive integer, then a
#' list of n rays.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayUniform <- function(n = NULL) {
  if (is.null(n)) {
    cartesianFromHorizontal(c(stats::runif(1, 0, 2 * pi), stats::runif(
      1,
      -1, 1
    )))
  } else {
    replicate(n, rayUniform(), simplify = FALSE)
  }
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineUniform <- function(n = NULL) {
  if (is.null(n)) {
    lower(rayUniform())
  } else {
    lapply(rayUniform(n), lower)
  }
}


#' Conversion from matrix to angle-axis representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians (between 0 and pi).
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' 
#' @noRd
#' 
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rotAntisymmetricFromVector <- function(v) {
  matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}

#' The ray that represents the line and that is closest to center.
#'
#' @param l A line (unit 3D vector).
#' @param center A ray (unit 3D vector).
#' 
#' @noRd
#' 
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineNearestRepresentative <- function(l, center) {
  if (dot(l, center) < 0) {
    -l
  } else {
    l
  }
}


#' Selecting an orientation representative near a given rotation.
#'
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A rotation matrix. The element of G R that is closest to center.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriNearestRepresentatives <- function(rs, center = rs[[1]], group) {
  lapply(rs, oriNearestRepresentative, center, group)
}




#' When trivial symmetry is used, the orientations are simply rotations. This case is so important that we have separate code, in rot.R, for doing it. But let's include the trivial group, for completeness.
#' @noRd
#' @keywords internal
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriTrivialGroup <- function() list(diag(c(1, 1, 1)))

#' Ray-in-plane symmetry is applicable to faults-with-slip-directions, such as slickensides (and certain minerals).
#' @noRd
#' @keywords internal
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriRayInPlaneGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    diag(c(-1, 1, -1))
  )
}
#' Line-in-plane symmetry is applicable to foliation-lineations, cylindrical fold orientations, triaxial ellipsoid orientations, and earthquake focal mechanisms (and olivine).
#' @noRd
#' @keywords internal
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriLineInPlaneGroup <- function() {
  list(
    diag(c(1, 1, 1)),
    diag(c(1, -1, -1)),
    diag(c(-1, 1, -1)),
    diag(c(-1, -1, 1))
  )
}
#' Trigonal trapezohedral is the point group of alpha-quartz.
#' @noRd
#' @keywords internal
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @noRd
#' @keywords internals
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rotMatrixAboutZ <- function(a) {
  cosine <- cos(a)
  sine <- sin(a)
  matrix(c(cosine, sine, 0, -sine, cosine, 0, 0, 0, 1), 3, 3)
}


rotMatrixFromAxisAngle <- function(ua) {
  m <- rotAntisymmetricFromVector(ua[1:3])
  diag(c(1, 1, 1)) + sin(ua[[4]]) * m + (1 - cos(ua[[4]])) * (m %*% m)
}

rotAntisymmetricFromVector <- function(v) {
  matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}




#' Sanitizing an orientation of a fault-with-slip.
#'
#' The issue here is that geologists like describing fault slips in terms of the
#' movement direction of the hanging wall, but the hanging wall is undefined for
#' vertical faults. A much better-behaved convention is to describe the
#' direction of slip using a vorticity vector. Inverse to
#' `geoPoleHangingFromPoleVorticity()`.
#'
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the
#' pole to a fault plane. The second row is the movement direction of the
#' hanging wall. (The third row is the cross product of the other two.)
#'
#' @return A real 3x3 matrix (special orthogonal). The first row is the
#' downward-pointing pole to the fault. The second row is the vorticity vector
#' of the slip. (The third row is the cross product of the other two.)
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
geoPoleVorticityFromPoleHanging <- function(rot) {
  if (rot[[1, 3]] < 0) {
    rbind(rot[1, ], -rot[3, ], rot[2, ])
  } else {
    rbind(-rot[1, ], rot[3, ], rot[2, ])
  }
}

#' Desanitizing an orientation of a fault-with-slip.
#'
#' Inverse to `geoPoleVorticityFromPoleHanging()`.
#'
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole
#' to a fault plane. The second row is the vorticity vector of the fault slip.
#' (The third row is the cross product of the other two.)
#'
#' @return A real 3x3 matrix (special orthogonal). The first row is the
#' downward-pointing pole to the fault. The second row is the movement
#' direction of the hanging wall. (The third row is the cross product of the
#' other two.)
#'
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
geoPoleHangingFromPoleVorticity <- function(rot) {
  if (rot[[1, 3]] < 0) {
    rbind(rot[1, ], rot[3, ], -rot[2, ])
  } else {
    rbind(-rot[1, ], -rot[3, ], -rot[2, ])
  }
}




# Mean and Variance -------------------------------------------------------------------------

#' @source `geologyGeometry` by Davis, J.R.
arithmeticMean <- function(xs, weights = NULL) {
  if (is.null(weights)) {
    Reduce("+", xs) / length(xs)
  } else {
    total <- Reduce("+", mapply(function(w, x, SIMPLIFY = TRUE) {
      w * x
    }, weights, xs))
    total / sum(weights)
  }
}


#' Extrinsic mean and scalar scatter of a set of rays.
#'
#' Scatter varies between 0 (tight concentration) and 1 (wide dispersion).
#' This scatter is denoted 1 - R-bar in Mardia and Jupp (2000), p. 163.
#' Arguably the preferred measure of scatter should be
#' 2 * (1 - R-bar) or 1 - R-bar^2, but this function implements neither of those.
#'
#' @param us A list of rays.
#' @return A list with elements `$mean` (ray) and $scatter (a real number
#' between 0 and 1).
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayMeanScatter <- function(us) {
  u <- arithmeticMean(us)
  r <- sqrt(dot(u, u))
  list(mean = (u / r), scatter = (1 - r))
}


#' Extrinsic mean of a set of rays.
#'
#' Convenience shortcut for `rayMeanScatter(us)$mean`.
#' @param us A list of rays.
#'
#' @return A ray.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayProjectedMean <- function(us) {
  rayNormalized(arithmeticMean(us))
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineProjectedMean <- function(us) {
  eig <- lineMeanScatter(us)
  eig$vectors[, 1]
}


#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayProjectedMean <- function(us) rayNormalized(arithmeticMean(us))


#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineMeanScatter <- function(us) {
  tMatrices <- lapply(us, function(u) outer(u, u))
  tMatrix <- arithmeticMean(tMatrices)
  eig <- eigen(tMatrix, symmetric = TRUE)
  eig
}




#' The Frechet (geodesic L^2) mean of a set of rays.
#'
#' An interative algorithm for computing the geodesic mean --- the ray that
#' minimizes the geodesic variance. The iterations continue until error squared
#' of epsilon is achieved or numSteps iterations have been used. Try multiple
#' seeds, to improve your chances of finding the global optimum.
#'
#' @param us A list of rays.
#' @param numSeeds A real number (positive integer). How many us to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations
#' to use.
#'
#' @return A list consisting of `$variance` (a real number), $mean (a ray),
#' `$error` (an integer) and `$minEigenvalue` (a real number). error should be
#' 0 and minEigenvalue should be positive. Otherwise there was some problem in
#' the optimization. If error is non-zero, then try increasing `numSteps`.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayMeanVariance <- function(us, numSeeds = 5, numSteps = 100) {
  f <- function(phiTheta) {
    rayVariance(us, cartesianFromSpherical(c(1, phiTheta)))
  }
  seeds <- sample(us, numSeeds)
  best <- list(variance = (pi^2))
  for (seed in seeds) {
    sol <- stats::optim(sphericalFromCartesian(seed), f,
      hessian = TRUE,
      control = list(maxit = numSteps)
    )
    if (sol$value < best$variance) {
      eigvals <- eigen(sol$hessian,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      best <- list(variance = sol$value, mean = cartesianFromSpherical(c(
        1,
        sol$par
      )), error = sol$convergence, minEigenvalue = min(eigvals))
    }
  }
  best
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineMeanVariance <- function(us, numSeeds = 5L, numSteps = 100L) {
  f <- function(phiTheta) {
    lineVariance(us, cartesianFromSpherical(c(1, phiTheta)))
  }
  seeds <- sample(us, numSeeds)
  best <- list(variance = (pi^2))
  for (seed in seeds) {
    sol <- stats::optim(sphericalFromCartesian(seed), f,
      hessian = TRUE,
      control = list(maxit = numSteps)
    )
    if (sol$value < best$variance) {
      eigvals <- eigen(sol$hessian,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      best <- list(variance = sol$value, mean = cartesianFromSpherical(c(
        1,
        sol$par
      )), error = sol$convergence, minEigenvalue = min(eigvals))
    }
  }
  best
}



#' L^2 variance of a set of rays about a point.
#'
#' I'm not sure about the weighting on this. If you change it, change regression
#' too.
#'
#' @param us A list of rays.
#' @param center A ray. Usually some kind of mean of the us.
#' @return A real number (in the interval `[0, pi^2 / 2]`).
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayVariance <- function(us, center) {
  sum(sapply(us, function(u) rayDistance(u, center)^2)) / (2 * length(us))
}

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineVariance <- function(us, center) {
  sum(sapply(us, function(u) lineDistance(u, center)^2)) / (2 * length(us))
}


#' Distance between two points on the unit sphere.
#'
#' @param u A ray.
#' @param v A ray.
#'
#' @return A real number. The angular difference between the two rays, in
#' radians. Always between 0 and pi, inclusive.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayDistance <- function(u, v) arcCos(dot(u, v))

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineDistance <- function(u, v) arcCos(abs(dot(u, v)))












#' The Frechet (geodesic L^2) variance of a set of orientations.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Often the mean of the rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative). The variance of the points G R about the point G center.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriVariance <- function(rs, center, group) {
  dists <- sapply(rs, oriDistance, center, group)
  sum(dists^2) / (2 * length(rs))
}



#' The distance between two orientations as points in SO(3) / G.
#'
#' @param r,q A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @returns A real number (in the interval `[0, pi]`). The distance from G R to G R.
#' @noRd
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
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
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
oriDiameter <- function(rs, group) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) oriDistance(rs[[i]], rs[[j]], group)))
  }
  max(sapply(2:(length(rs)), f))
}


# Tangent space methods, e.g. PCA --------------------------------------------------------------------------

#' Principal component analysis in the tangent space.
#'
#' Appropriate only if the data set is tightly concentrated about the center.
#'
#' @param us A list of rays.
#' @param center A ray. Typically the geodesic mean of the us.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to
#' return on each of the four geodesics through the center.
#'
#' @return A list consisting of `$rotation` (3x3 rotation matrix),
#' `$magnitudes` (2D real vector, non-negative),
#' `$directions` (2x2 real matrix, whose columns are unit-length vectors),
#' `$pcsFromRay` (an R function from rays to 2-dimensional vectors), and
#' `$rayFromPCs` (an R function from 2-dimensional vectors to rays). The
#' `$magnitudes` are analogous to sample standard deviations. They are in
#' decreasing order. The `$directions` are the corresponding directions,
#' suitable for use in `rayPointFromTangentVector()` along with $rotation.
#' If numPoints >= 1, then there is also a `$curves` field (list of two lists
#' of (2 numPoints + 1) rays).
#'
#' @noRd
#'
#' @source `geologyGeometry` by Davis, J.R.
rayPrincipalComponentAnalysis <- function(us, center, numPoints = 0) {
  # Map the rays into the tangent space at the mean.
  perp <- rayOrthogonalUniform(center)
  rot <- rbind(center, perp, cross(center, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  # Compute the usual covariance stuff.
  covar <- arithmeticMean(lapply(vs, function(v) {
    outer(v, v)
  }))
  eig <- eigen(covar, symmetric = TRUE)
  mags <- sqrt(eig$values)
  dirs <- eig$vectors
  # Compute functions for transferring back and forth.
  pcsFromRay <- function(u) {
    as.numeric(t(dirs) %*% rayTangentVectorFromPoint(u, rot))
  }
  rayFromPCs <- function(w) {
    rayPointFromTangentVector(as.numeric(dirs %*% w), rot)
  }
  result <- list(
    rotation = rot, magnitudes = mags, directions = dirs, pcsFromRay = pcsFromRay,
    rayFromPCs = rayFromPCs
  )
  # Include points for visualization, if desired.
  if (numPoints >= 1) {
    f <- function(s, magDir) {
      rayPointFromTangentVector(s / numPoints * magDir, rot)
    }
    curve1 <- lapply(-numPoints:numPoints, f, mags[[1]] * dirs[, 1])
    curve2 <- lapply(-numPoints:numPoints, f, mags[[2]] * dirs[, 2])
    result$curves <- list(curve1, curve2)
  }
  result
}

linePrincipalComponentAnalysis <- function(us, center, numPoints = 0) {
  rs <- lapply(us, lineNearestRepresentative, center)
  pca <- rayPrincipalComponentAnalysis(rs, center, numPoints)
  result <- list(
    rotation = pca$rotation, magnitudes = pca$magnitudes,
    directions = pca$directions, pcsFromLine = pca$pcsFromRay,
    lineFromPCs = pca$rayFromPCs
  )
  if (numPoints >= 1) {
    result$curves <- pca$curves
  }
  result
}


#' Principal geodesic analysis in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to return on each of the six geodesics through the center.
#' @return A list consisting of $magnitudes (3D real vector, nonnegative) and $directions (3x3 real matrix, whose columns are unit-length vectors). The $magnitudes are analogous to sample standard deviations. They are in decreasing order. The $directions are the corresponding directions, suitable for use in rotMatrixFromLeftTangent. If numPoints >= 1, then there is also a $curves field (list of three lists of (2 numPoints + 1) rotation matrices).
#' @noRd
#' @source `geologyGeometry` by Davis, J.R.
rotLeftPrincipalComponentAnalysis <- function(rs, center, numPoints = 0) {
  eig <- eigen(rotLeftCovariance(rs, center), symmetric = TRUE)
  mags <- sqrt(eig$values)
  dirs <- eig$vectors
  result <- list(magnitudes = mags, directions = dirs)
  if (numPoints >= 1) {
    f <- function(s, magDir) {
      rotMatrixFromLeftTangent(s / numPoints * magDir, center)
    }
    curve1 <- lapply(-numPoints:numPoints, f, mags[[1]] * dirs[, 1])
    curve2 <- lapply(-numPoints:numPoints, f, mags[[2]] * dirs[, 2])
    curve3 <- lapply(-numPoints:numPoints, f, mags[[3]] * dirs[, 3])
    result$curves <- list(curve1, curve2, curve3)
  }
  result
}


#' Principal geodesic analysis in the tangent space at a given orientation.
#'
#' Appropriate only if the sample is tightly concentrated near the center. In 
#' fact, this function is a thin wrapper for rotLeftPrincipalComponentAnalysis.
#' 
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the geodesic mean of the rs.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to 
#' return on each of the six geodesics through the center.
#' 
#' @return A list consisting of $magnitudes (3D real vector, nonnegative) and 
#' $directions (3x3 real matrix, whose columns are unit-length vectors). The 
#' $magnitudes are analogous to sample standard deviations. They are in 
#' decreasing order. The $directions are the corresponding directions, suitable 
#' for use in rotMatrixFromLeftTangent. If numPoints >= 1, then there is also a 
#' $curves field (list of three lists of (2 numPoints + 1) rotation matrices).
#' 
#' @noRd
#' 
#' @source `geologyGeometry` by Davis, J.R.
oriLeftPrincipalComponentAnalysis <- function(rs, center, group, numPoints = 0) {
  rs <- lapply(rs, oriNearestRepresentative, center, group)
  rotLeftPrincipalComponentAnalysis(rs, center, numPoints)
}



#' Unwrapping unit sphere into a tangent plane.
#'
#' Inverse to `rayPointFromTangentVector()`.
#'
#' @param q A ray.
#' @param rotation A 3x3 real matrix (special orthogonal).
#'
#' @return A 2-dimensional real vector.
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @noRd
rayTangentVectorFromPoint <- function(q, rotation) {
  v <- rayLog(rotation[1, ], q)
  w <- as.numeric(rotation %*% v)[2:3]
  w
  # Does this faster code do the same thing?
  # as.numeric(rotation %*% q)[2:3]
}

#' Wrapping a plane around the unit sphere.
#'
#' This function composes the exponential map with a non-canonical isomorphism
#' to the plane \eqn{R^2}. The isomorphism is specified by the user through a
#' rotation matrix `R`. The first row of `R` is regarded as a point p on the unit
#' sphere. The other two rows define an isomorphism between \eqn{R^2} and the tangent
#' plane to the unit sphere at p. `w` is mapped through this isomorphism into the
#' tangent plane, and then into the sphere via the exponential map. Inverse to
#' `rayTangentVectorFromPoint()`.
#'
#' @param w A 2-dimensional real vector.
#' @param rotation A 3x3 real matrix (special orthogonal).
#'
#' @return A ray.
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @noRd
rayPointFromTangentVector <- function(w, rotation) {
  v <- as.numeric(t(rotation) %*% c(0, w))
  q <- rayExp(rotation[1, ], v)
  q
  # Does this other code do the same thing?
  # v <- c(sqrt(1 - vec[[1]]^2 - vec[[2]]^2), vec)
  # as.numeric(t(rotation) %*% v)
}



#' Sample covariance matrix, approximated in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the geodesic mean of the `rs`.
#' 
#' @return A 3x3 real matrix (symmetric, non-negative eigenvalues).
#' 
#' @export
rotLeftCovariance <- function(rs, center) {
  vs <- lapply(rs, rotLeftTangentFromMatrix, center)
  ms <- lapply(vs, function(v) {
    outer(v, v)
  })
  arithmeticMean(ms)
}

#' Maps the space of rotations into one of its tangent spaces.
#'
#' Converts rotation matrix \eqn{R = \exp(V) C} into vector v. Appropriate only 
#' if the two given rotations are close to each other.
#' 
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' 
#' @return A 3D real vector.
#' 
#' @noRd
rotLeftTangentFromMatrix <- function(r, center) {
  ua <- rotAxisAngleFromMatrix(t(center) %*% r)
  ua[1:3] * ua[[4]]
}

#' Conversion from matrix to angle-axis representation.
#'
#' @param r A rotation matrix.
#' 
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle 
#' in radians (between 0 and pi).
#' 
#' @noRd
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

#' Maps a tangent space into the space of rotations.
#'
#' Converts tangent vector `v` into rotation matrix \eqn{R = C \exp(V)}.
#' 
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' 
#' @return A 3x3 rotation matrix.
#' @noRd
rotMatrixFromLeftTangent <- function(v, center) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0) {
    center
  } else {
    center %*% rotMatrixFromAxisAngle(c(v / nrm, nrm))
  }
}


# Inference --------------------------------------------------------------------

#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
lineBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls,
    length(ls),
    replace = TRUE
  )), simplify = FALSE)
  bootMean <- lineProjectedMean(boots)
  boots <- lapply(boots, function(u) {
    if (dot(u, bootMean) >= 0) {
      u
    } else {
      -u
    }
  })
  inf <- rayMahalanobisPercentiles(boots, center = bootMean, ...)

  inf$pvalueLine <- function(l) {
    if (is.spherical(l)) l <- unclass(Vec3(l))
    if (dot(l, inf$center) < 0) {
      inf$pvalue(-l)
    } else {
      inf$pvalue(l)
    }
  }
  inf
}

#' Bootstrapped extrinsic mean with percentile confidence region and hypothesis
#' tests.
#'
#' The inference is based on percentiles of Mahalanobis distance in the tangent
#' space at the mean of the bootstrapped means. The user should check that the
#' bootstrapped means form a tight ellipsoidal cluster, before taking such a
#' region seriously.
#'
#' @param ls A list of rays.
#' @param numBoots A real number (positive integer). The number of
#' bootstrapped means to compute. 10,000 might be a good value.
#' @param ... Other arguments to be passed to the underlying
#' `rayMahalanobisPercentiles()` function.
#'
#' @return A list. See `rayMahalanobisPercentiles()`.
#'
#' @noRd
#'
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
rayBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls,
    length(ls),
    replace = TRUE
  )), simplify = FALSE)
  bootMean <- rayProjectedMean(boots)
  rayMahalanobisPercentiles(boots, bootMean, ...)
}

#' Elliptical region based on percentiles of Mahalanobis distance in the tangent
#' space.
#'
#' @param us A list of rays.
#' @param center A ray. Usually something like the geodesic mean of us.
#' @param alpha A real number, which is assumed to be 0.05 unless you specify 0.01.
#' @param numPoints A real number (non-negative integer). The resolution with
#' which to sample the boundary curve of the (1 - alpha) * 100% percentile region.
#' @param doIsotropic Logical. If TRUE, forces the inverse covariance to the
#' identity matrix and hence the region to be circular.
#'
#' @return A list. `$us` is `us`. `$center` is center. `$covarInv` is the inverse
#' covariance matrix in the tangent space, which is just the identity if
#' `doIsotropic` is `TRUE`. rotation is the rotation used in
#' `rayTangentVectorFromPoint()`, etc.
#' `$q000`, `$q025`, `$q050`, `$q075`, `$q095`, `$q099` are quantiles of
#' Mahalanobis norm.
#' `$pvalue` is an R function that assigns to any ray its p-value, meaning the
#' fraction of us that are farther from center than the given ray.
#' `$alpha` is alpha. $angles is a pair of two angles, in radians, giving the
#' semi-axis lengths of the confidence ellipse. If `numPoints > 0`, then the
#' list also has an element $points.
#' `$points` is a list of `numPoints + 1` rays describing the boundary of the
#' region, in order, with the first and last points identical.
#'
#' @noRd
#'
#' @source `geologyGeometry` by Davis, J.R.
rayMahalanobisPercentiles <- function(us, center, alpha = 0.05, numPoints = 0, doIsotropic = FALSE) {
  perp <- rayOrthogonalUniform(center)
  rot <- rbind(center, perp, cross(center, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  if (doIsotropic) {
    covarInv <- diag(c(1, 1))
  } else {
    covar <- arithmeticMean(lapply(vs, function(v) {
      outer(v, v)
    }))
    covarInv <- solve(covar)
  }
  norms <- sapply(vs, function(v) {
    sqrt(v %*% covarInv %*% v)
  })
  empiricalCDF <- stats::ecdf(norms)
  f <- function(u) {
    if (is.spherical(u)) u <- unclass(Vec3(u))
    v <- rayTangentVectorFromPoint(u, rot)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  qs <- stats::quantile(norms, probs = c(
    0, 0.25, 0.5, 0.75, 0.95,
    0.99, 1
  ), names = FALSE)
  if (alpha != 0.05 && alpha != 0.01) {
    alpha <- 0.05
  }
  if (alpha == 0.01) {
    q <- qs[[6]]
  } else {
    q <- qs[[5]]
  }
  eig <- eigen(covarInv, symmetric = TRUE)
  radii <- q * eig$values^(-1 / 2)
  result <- list(
    us = us, pvalue = f, center = center, covarInv = covarInv,
    rotation = rot, q000 = qs[[1]], q025 = qs[[2]], q050 = qs[[3]],
    q075 = qs[[4]], q095 = qs[[5]], q099 = qs[[6]], q100 = qs[[7]],
    alpha = alpha, angles = radii
  )
  if (numPoints > 0) {
    circle <- lapply(0:numPoints, function(s) {
      c(cos(s *
        2 * pi / numPoints), sin(s * 2 * pi / numPoints))
    })
    vs <- lapply(circle, function(v) {
      as.numeric(eig$vectors %*%
        (radii * v))
    })
    result$points <- lapply(
      vs, rayPointFromTangentVector,
      rot
    )
    result$alpha <- alpha
  }
  result
}







# Regression -------------------------------------------------------------------

rayRegressionSmallCircle <- function(xs, us, numSeeds = 5, numSteps = 1000, numPoints = 0, angleBound = Inf) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      r <- rotMatrixFromAxisAngle(c(pole, x * phiThetaAlphaTauSigma[[3]]))
      as.numeric(r %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(rayDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value = (2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(
      sphericalFromCartesian(rayUniform())[2:3],
      stats::runif(1, -pi, pi),
      sphericalFromCartesian(rayUniform())[2:3]
    )
    if (angleBound == Inf) {
      solution <- stats::optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B"
      )
    } else {
      solution <- stats::optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B",
        lower = c(-Inf, -Inf, -angleBound, -Inf, -Inf),
        upper = c(Inf, Inf, angleBound, Inf, Inf)
      )
    }
    if (solution$value <= best$value) {
      best <- solution
    }
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric = TRUE, only.values = TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- rayVariance(us, rayProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {
    as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  }
  results <- list(
    error = best$convergence, minEigenvalue = min(eigvals), pole = pole, angle = angle,
    rSquared = rSquared, prediction = pred
  )
  if (numPoints >= 1) {
    ys <- seq(from = min(xs), to = max(xs), length.out = numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

# Warning: Not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
rayRegressionSmallCircleRescaled <- function(xs, ls, angleBound = Inf, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- rayRegressionSmallCircle(
    scales(xs), ls,
    angleBound = (angleBound * (x1 - x0)), ...
  )
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {
    results$rescaledPrediction((x - x0) / (x1 - x0))
  }
  results$angle <- results$angle / (x1 - x0)
  results
}

lineRegressionSmallCircle <- function(xs, us, numSeeds = 5, numSteps = 1000, numPoints = 0, angleBound = Inf) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      angle <- x * phiThetaAlphaTauSigma[[3]]
      as.numeric(rotMatrixFromAxisAngle(c(pole, angle)) %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(lineDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value = (2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(
      sphericalFromCartesian(lineUniform())[2:3],
      stats::runif(1, -pi, pi),
      sphericalFromCartesian(lineUniform())[2:3]
    )
    if (angleBound == Inf) {
      solution <- stats::optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B"
      )
    } else {
      solution <- stats::optim(
        seed, f,
        hessian = TRUE, control = list(maxit = numSteps), method = "L-BFGS-B",
        lower = c(-Inf, -Inf, -angleBound, -Inf, -Inf),
        upper = c(Inf, Inf, angleBound, Inf, Inf)
      )
    }
    if (solution$value <= best$value) {
      best <- solution
    }
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric = TRUE, only.values = TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- lineVariance(us, lineProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {
    as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  }
  results <- list(
    error = best$convergence, minEigenvalue = min(eigvals), pole = pole, angle = angle,
    rSquared = rSquared, prediction = pred
  )
  if (numPoints >= 1) {
    ys <- seq(from = min(xs), to = max(xs), length.out = numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

# Not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
lineRegressionSmallCircleRescaled <- function(xs, us, angleBound = Inf, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- lineRegressionSmallCircle(
    xs = scales(xs), us = us, angleBound = (angleBound * (x1 - x0)), ...
  )
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {
    results$rescaledPrediction((x - x0) / (x1 - x0))
  }
  results$angle <- results$angle / (x1 - x0)
  results
}
