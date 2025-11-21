#' The Fr&#233;chet (geodesic \eqn{L^2}) variance
#'
#' Dispersion  measured using the Fr&#233;chet variance, i.e the sum of the squared
#' geodesic distances between all vectors and a specified vector.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"`,
#' where the rows are the observations and the columns are the coordinates.
#' @param y Only for variance. object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"` about which the Fr&#233;chet variance should be calculated for.
#' If `NULL` (the default), Fr&#233;chet variance about the Fr&#233;chet mean.
#' @param ... parameters passed to [geodesic_meanvariance_ray()] (if `x` is a Ray), [geodesic_meanvariance_line()] (if `x` is a Vec3, Line or Plane)
#' or [geodesic_mean_pair()] (if `x` is a Pair or a Fault).
#' @inheritParams geodesic_mean_pair
#'
#' @details The variance of a dataset \eqn{{x_1, \ldots, x_n}} about a vector \eqn{y} is defined as
#' \deqn{ \Psi(x) = \frac{1}{2n} \sum_{i=1}^n d_G(y, x_i)^2}
#' where \eqn{d_G(x, y)} is the geodesic distance between vectors \eqn{x} and \eqn{y} (see [angle()]).
#'
#' @name geodesic-var
#'
#' @returns the Fr&#233;chet variance as a numeric number. Because distances in SO(3) never exceed \eqn{\pi}, the maximum possible variance
#' is \eqn{\frac{\pi^2}{2} \approx 4.93}.
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. \doi{10.1016/j.jsg.2017.01.002}
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @seealso [geodesic_mean()] for the Fr&#233;chet mean, [sph_mean()] for the arithmetic mean, [projected_mean()] for projected mean
#'
#' @examples
#' set.seed(20250411)
#' geodesic_var(example_planes, example_planes[1, ])
#' geodesic_var(example_planes)
NULL

#' @rdname geodesic-var
#' @export
geodesic_var <- function(x, ...) UseMethod("geodesic_var")

#' @rdname geodesic-var
#' @export
geodesic_var.Vec3 <- function(x, y = NULL, ...) {
  if (is.null(y)) {
    geodesic_var_line(x, ...)
  } else {
    stopifnot(is.spherical(y) & nrow(y) == 1)
    lineVariance(vec_list(x), as.vector(unclass(Vec3(y))))
  }
}

#' @rdname geodesic-var
#' @export
geodesic_var.Line <- function(x, y = NULL, ...) {
  if (is.null(y)) {
    geodesic_var_line(x, ...)
  } else {
    stopifnot(is.spherical(y) & nrow(y) == 1)
    lineVariance(vec_list(x), as.vector(unclass(Vec3(y))))
  }
}

#' @rdname geodesic-var
#' @export
geodesic_var.Plane <- function(x, y = NULL, ...) {
  if (is.null(y)) {
    geodesic_var_line(x, ...)
  } else {
    stopifnot(is.spherical(y) & nrow(y) == 1)
    lineVariance(vec_list(x), as.vector(unclass(Vec3(y))))
  }
}

#' @rdname geodesic-var
#' @export
geodesic_var.Ray <- function(x, y = NULL, ...) {
  if (is.null(y)) {
    geodesic_var_ray(x, ...)
  } else {
    stopifnot(is.spherical(y) & nrow(y) == 1)
    rayVariance(vec_list(x), as.vector(unclass(Vec3(y))))
  }
}

#' @rdname geodesic-var
#' @export
geodesic_var.Pair <- function(x, y = NULL, group = NULL, ...) {
  if (is.null(y)) {
    geodesic_var_pair(x, ...)
  } else {
    stopifnot(is.spherical(y) & nrow(y) == 1)

    if (is.null(group)) {
      group <- if (inherits(y, "Fault")) "triclinic" else "triclinic"
    }

    group_mat <- symmetry_group(group)

    oriVariance(vec_list(x), pair2rot(y), group = group_mat)
  }
}


#' The Fr&#233;chet (geodesic \eqn{L^2}) mean
#'
#' An iterative algorithm for computing the Fr&#233;chet mean, i.e. the vector that
#' minimizes the Fr&#233;chet variance.
#'
#' @inheritParams geodesic-var
#'
#' @name geodesic-mean
#'
#' @returns the Fr&#233;chet mean vector as an object of class `x`.
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. \doi{10.1016/j.jsg.2017.01.002}
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#'
#' @seealso [geodesic_var()] for Fr&#233;chet variance, [sph_mean()] for the arithmetic mean, [projected_mean()] for projected mean
#' @source `geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/
#' @examples
#' set.seed(20250411)
#' geodesic_mean(example_planes)
NULL

#' @rdname geodesic-mean
#' @export
geodesic_mean <- function(x, ...) UseMethod("geodesic_mean")

#' @rdname geodesic-mean
#' @export
geodesic_mean.Vec3 <- function(x, ...) geodesic_mean_line(x, ...)

#' @rdname geodesic-mean
#' @export
geodesic_mean.Line <- function(x, ...) geodesic_mean_line(x, ...)

#' @rdname geodesic-mean
#' @export
geodesic_mean.Ray <- function(x, ...) geodesic_mean_ray(x, ...)

#' @rdname geodesic-mean
#' @export
geodesic_mean.Plane <- function(x, ...) geodesic_mean_line(x, ...)

#' @rdname geodesic-mean
#' @export
geodesic_mean.Pair <- function(x, ...) geodesic_mean_pair(x, ...)










#' The Fr&#233;chet (geodesic \eqn{L^2}) mean of a set of lines or rays
#'
#' An iterative algorithm for computing the Fr&#233;chet mean — the line or ray that
#' minimizes the Fr&#233;chet variance. The iterations continue until error squared of
#' epsilon is achieved or `steps` iterations have been used. Try multiple
#' seeds, to improve your chances of finding the global optimum.
#'
#' @inheritParams sph_mean
#' @param seeds positive integer. How many `x` to try as seeds
#' @param steps positive integer. Bound on how many iterations to use.
#' @param ... parameters passed to [geodesic_meanvariance_line()] or [geodesic_meanvariance_ray()]
#'
#' @returns `geodesic_meanvariance_*` (* either line or ray) returns a `list` consisting of
#' `$variance` (numeric), `$mean` (a line),
#' `$error` (an integer) and `$min.eigenvalue` (numeric).
#' `geodesic_mean_*` and `geodesic_var_*` are convenience wrapper and only
#' return the mean and the variance, respectively.
#'
#' @details Error should be `0` and `min.eigenvalue` should be positive.
#' Otherwise there was some problem in the optimization. If error is non-zero, then try increasing `steps`.
#' @name geodesic-line
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. https://doi.org/10.1016/j.jsg.2017.01.002
#' @source lineMeanVariance from geologyGeometry (J.R. Davis)
#' @seealso [sph_mean()] for arithmetic mean and [geodesic_mean_pair()] for geodesic mean of pairs.
#'
#' @examples
#' geodesic_meanvariance_line(example_lines)
#' geodesic_mean_line(example_lines)
#' geodesic_var_line(example_lines)
NULL

#' @rdname geodesic-line
#' @export
geodesic_mean_line <- function(x, ...) geodesic_meanvariance_line(x, ...)$mean

#' @rdname geodesic-line
#' @export
geodesic_mean_ray <- function(x, ...) geodesic_meanvariance_ray(x, ...)$mean


#' @rdname geodesic-line
#' @export
geodesic_var_line <- function(x, ...) geodesic_meanvariance_line(x, ...)$variance

#' @rdname geodesic-line
#' @export
geodesic_var_ray <- function(x, ...) geodesic_meanvariance_ray(x, ...)$variance

#' @rdname geodesic-line
#' @export
geodesic_meanvariance_line <- function(x, seeds = 5L, steps = 100L) {
  stopifnot(is.Vec3(x) | is.Line(x) | is.Plane(x))
  res <- lineMeanVariance(vec_list(x), numSeeds = seeds, numSteps = steps)
  m <- res$mean |> as.Vec3()

  names(res) <- c("variance", "mean", "error", "min.eigenvalue")

  res$mean <- Spherical(m, class(x)[1])
  return(res)
}

#' @rdname geodesic-line
#' @export
geodesic_meanvariance_ray <- function(x, seeds = 5L, steps = 100L) {
  stopifnot(is.Vec3(x) | is.Ray(x) | is.Plane(x))
  res <- rayMeanVariance(vec_list(x), numSeeds = seeds, numSteps = steps)
  m <- res$mean |> as.Vec3()

  names(res) <- c("variance", "mean", "error", "min.eigenvalue")

  res$mean <- Spherical(m, class(x)[1])
  return(res)
}







#' Mean orientation of a set of pairs or faults
#'
#' The Frechet (geodesic \eqn{L^2}) mean and variance of a pair of foliations
#' and lineations
#'
#' @param x object of class `"Pair"` or `"Fault"`
#' @param group Symmetry group of `x`. See [symmetry_group()] for details
#' If `NULL`, the group will be
#' automatically picked based on the class of `x`.
#'
#' @returns object of class `"Pair"` or `"Fault"`, respectively
#' @name mean-pair
#' @seealso [sph_mean()] for arithmetic mean, and [geodesic_mean_line()] for geodesic mean of lines.
#'
#' @references Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis
#' for three-dimensional orientational data. Journal of Structural Geology,
#' 96, 65–89. https://doi.org/10.1016/j.jsg.2017.01.002
#'
#' @source oriMeanVariance from geologyGeometry (J.R. Davis)
#'
#' @examples
#' my_fault <- Fault(
#'   c("a" = 120, "b" = 120, "c" = 100),
#'   c(60, 60, 50),
#'   c(110, 25, 30),
#'   c(58, 9, 23),
#'   c(1, -1, 1)
#' )
#' geodesic_mean_pair(my_fault)
#' geodesic_var_pair(my_fault)
NULL

#' @rdname mean-pair
#' @export
geodesic_mean_pair <- function(x, group = NULL) {
  res <- .sph_meanvar_pair(x, group)
  rot2pair(res$mean, fault = inherits(x, "Fault"))
}

#' @rdname mean-pair
#' @export
geodesic_var_pair <- function(x, group = NULL) {
  .sph_meanvar_pair(x, group)$variance
}

.sph_meanvar_pair <- function(x, group) {
  xvec <- pair2rot(x)

  if (is.null(group)) {
    group <- if (inherits(x, "Fault")) "triclinic" else "monoclinic"
  }

  group_mat <- symmetry_group(group)
  ori_mean_variance(xvec, group = group_mat)
}


#' Orientation matrix from fault planes and slip directions
#'
#' Converts a set of planes and lines into a list of rotation matrices.
#'
#' @param x object of class `"Pair"` or `"Fault"`
#' @noRd
#' @returns list of rotation matrices
#'
#' @examples
#' my_fault <- Fault(
#'   c("a" = 120, "b" = 120, "c" = 100),
#'   c(60, 60, 50),
#'   c(110, 25, 30),
#'   c(58, 9, 23),
#'   c(1, -1, 1)
#' )
#' pair2rot(my_fault)
pair2rot <- function(x) {
  stopifnot(is.Pair(x))
  p <- Fault_plane(x) |> Vec3()
  l <- Fault_slip(x) |> Vec3()

  if (inherits(x, "Fault")) {
    l <- x[, 5] * l
  }
  cross <- crossprod.Vec3(p, l)

  rotm <- rot_projected_matrix(list(pole = p, direction = l, cross = cross)) |>
    as.Rotation()
  return(rotm)
}

is.Rotation <- function(x) inherits(x, "Rotation")

as.Rotation <- function(x) {
  class(x) <- append("Rotation", class(x))
  return(x)
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
#' my_fault <- Fault(
#'   c("a" = 120, "b" = 120, "c" = 100),
#'   c(60, 60, 50),
#'   c(110, 25, 30),
#'   c(58, 9, 23),
#'   c(1, -1, 1)
#' )
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





### SYMMETRY GROUPS ###

#' Symmetry groups
#'
#' Frequently orientations are subject to some symmetry group \eqn{\mathbb{G}}, which is a
#' finite set of rotations (satisfying certain properties). For any rotation \eqn{Q}
#' in \eqn{\mathbb{G}}, the rotation matrix \eqn{Q \cdot R} represents the same orientation as \eqn{R} does.
#'
#' @param group character (symmetry class) or integer number of the enantiomorphic point group. See table below for details.
#'  Also accepts `""ray_in_plane"` (equivalent to triclinic symmetry) and `"line_in_plane"` (monoclinic).
#'
#' @details
#' \tabular{lcl}{
#' Symmetry \tab Enantiomorphic Point Group \tab Example \cr
#' `triclinic` 	\tab 1 \tab Ray in plane, plagioclase \cr
#' `monoclinic` \tab 2 \tab Line in plane, orthoclase, gypsum, muscovite, clinopyroxene, clinoamphibole \cr
#' `orthorhombic` \tab 222 \tab olivine, aragonite, marcasite, orthopyroxenes \cr
#' `tetragonal` \tab 	4 \tab Pyramidal: zircon \cr
#'            \tab 422 \tab Trapezohedral\cr
#' `trigonal` \tab 	3 \tab Pyramidal, Rhombohedral \cr
#'          \tab	32 \tab Trapezohedral: \eqn{\alpha}-Quartz \cr
#' `hexagonal` \tab	6 \tab Pyramidal, Rhombohedral\cr
#'           \tab	622 \tab Trapezohedral: \eqn{\beta}-Quartz \cr
#' `cubic` \tab	 23 \tab Tetartoidal \cr
#'       \tab	432 \tab Hexoctahedral: galena, pyrite, fluorite \cr
#' }
#'
#' @return list of rotation parameters
#'
#' @export
#'
#' @examples
#' symmetry_group("triclinic")
#'
#' symmetry_group(2)
symmetry_group <- function(group = c("triclinic", "ray_in_plane", "line_in_plane", "monoclinic", "trigonal-trapezohedral", "hexagonal-trapezohedral", "trivial")) {
  if (is.character(group)) {
    group <- match.arg(group)

    switch(group,
      "ray_in_plane" = oriRayInPlaneGroup(), # 1 triclinic
      "triclinic" = oriRayInPlaneGroup(),
      "line_in_plane" = oriLineInPlaneGroup(), # 2 monoclinic
      "monoclinie" = oriLineInPlaneGroup(), # 2 monoclinic
      "orthorhombic" = NULL,
      "tetragonal" = NULL,
      "trigonal-trapezohedral" = oriTrigonalTrapezohedralGroup(), # 32 trigonal-trapezohedral
      "hexagonal-trapezohedral" = oriHexagonalTrapezohedralGroup(), # 622 hexagonal-trapezohedral
      "cubic" = NULL,
      "trivial" = oriTrivialGroup()
    )
  } else {
    stopifnot(group %in% c(1, 2, 222, 4, 422, 3, 32, 6, 622, 23, 432))
    switch(as.character(group),
      "1" = oriRayInPlaneGroup(), # 1 triclinic
      "2" = oriLineInPlaneGroup(), # 2 monoclinic
      "222" = NULL,
      "4" = NULL,
      "422" = NULL,
      "3" = NULL,
      "32" = oriTrigonalTrapezohedralGroup(), # 32 trigonal-trapezohedral
      "6" = NULL,
      "622" = oriHexagonalTrapezohedralGroup(), # 622 hexagonal-trapezohedral
      "23" = NULL,
      "432" = NULL
    )
  }
}
