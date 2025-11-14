# PCA --------------------------------------------------------------------------


#' Principal component (geodesic) analysis in the tangent space.
#' 
#' @inheritParams geodesic-mean
#' @param center A spherical object. Typically the geodesic mean of the x.
#' @param n  real number (integer, 0 or >= 3). The number of points to
#' return on each of the four geodesics through the center.
#' @inheritParams geodesic_mean_pair
#'
#' @details Appropriate only if the data set is tightly concentrated about the center.
#'
#' @return A list consisting of
#' \describe{
#' \item{`rotation`}{3x3 rotation matrix}
#' \item{`magnitudes`}{2D real vector, non-negative. Magnitudes are
#' analogous to sample standard deviations. They are in  decreasing order and
#' quantify how dispersed the data are in those two directions.}
#' \item{`directions`}{2x2 real matrix, whose columns are unit-length vectors.
#' The corresponding directions to the magnitudes}
#' \item{`pcsFromRay`}{function to convert rays to 2-dimensional vectors}
#' \item{`rayFromPCs`}{function to convert 2-dimensional vectors to rays}
#' \item{`curves`}{list of two lists of (2 `n` + 1) rays, only if `n` >= 1}
#' \item{`tangents`}{Tangents from `directions` and `rotation`}
#' }
#' 
#' If `x` is a `"Pair"` object the list only contains `magnitudes`, `directions` 
#' and `curves` (if `n`>=1).
#'
#' @name sph-prcomp
#'
#' @examples
#' res <- prcomp(example_lines, n = 10)
#'
#' stereoplot(sub = paste("SD1:", round(res$magnitudes[1], 2), "| SD2:",round(res$magnitudes[2], 2)))
#' points(example_lines, pch = 16, cex = .8)
#' invisible(lapply(res$curves, stereo_lines, col = "red", lwd = 1.5))
NULL

#' @rdname sph-prcomp
#' @exportS3Method stats::prcomp
prcomp.Vec3 <- function(x, center = geodesic_mean(x), n = 0L) {
  center <- Vec3(center) |>
    unclass() |>
    as.vector()

  res <- vec_list(x) |>
    rayPrincipalComponentAnalysis(center, numPoints = n)

  res$tangents <- rbind(
    rayPointFromTangentVector(res$directions[, 1], res$rotation),
    rayPointFromTangentVector(res$directions[, 2], res$rotation)
  ) |> as.Vec3()

  res$rotation <- as.Rotation(res$rotation )

  # res$direction <- list_vec(res$direction)
  if (n > 0) {
    res$curves <- lapply(res$curves, function(x) {
      do.call(rbind, x) |>
        as.Vec3()
    })
  }

  return(res)
}

#' @rdname sph-prcomp
#' @exportS3Method stats::prcomp
prcomp.Ray <- function(x, center = geodesic_mean(x), n = 0L) {
  res <- pcomp.Vec3(x)

  # res$direction <- Ray(res$direction)
  # if(n > 0) res$curves <- lapply(res$curves, Ray)

  res$tangents <- Ray(res$tangents)

  return(res)
}

#' @rdname sph-prcomp
#' @exportS3Method stats::prcomp
prcomp.Line <- function(x, center = geodesic_mean(x), n = 0L) {
  center <- Vec3(center) |>
    unclass() |>
    as.vector()

  res <- vec_list(x) |>
    linePrincipalComponentAnalysis(center, numPoints = n)

  res$tangents <- rbind(
    rayPointFromTangentVector(res$directions[, 1], res$rotation),
    rayPointFromTangentVector(res$directions[, 2], res$rotation)
  ) |>
    as.Vec3() |>
    Line()

  res$rotation <- as.Rotation(res$rotation)

  # res$direction <- list_vec(res$direction) |> Line()
  if (n > 0) {
    res$curves <- lapply(res$curves, function(x) {
      do.call(rbind, x) |>
        as.Vec3()
    })
  }

  return(res)
}

#' @rdname sph-prcomp
#' @exportS3Method stats::prcomp
prcomp.Plane <- function(x, center = geodesic_mean(x), n = 0L) {
  res <- pcomp.Line(x)

  # res$direction <- Ray(res$direction)
  # if(n > 0) res$curves <- lapply(res$curves, Ray)

  res$tangents <- Plane(res$tangents)

  return(res)
}

#' @rdname sph-prcomp
#' @exportS3Method stats::prcomp
prcomp.Pair <- function(x, center = geodesic_mean(x), n = 0L, group = NULL) {
  xl <- pair2rot(x)
  center <- pair2rot(center)[[1]] |> unclass()

  res <- if(is.null(group)) {
    rotLeftPrincipalComponentAnalysis(xl, center, numPoints = n)
  } else {
    oriLeftPrincipalComponentAnalysis(xl, center, numPoints = n, group = group)
  }

  if (n > 0) {
    res$curves <- lapply(res$curves, function(i) {
      lapply(i, function(j) {
        as.Rotation(j) |>
          rot2pair() |>
          unclass()
      })
    }) |>
      lapply(function(i) {
        as.Pair(do.call(rbind, i))
      })
  }

  return(res)
}
