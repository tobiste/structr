#' Plot spherical densities in a stereonet
#'
#' Kamb counts and densities on the sphere. `contour` plots the contour lines,
#' `contourf` displays a contour plot with the areas between the contours filled,
#' and `image` creates a dense grid of colored rectangles.
#'
#' @param x object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"` or `'sph_density'` (for plotting only).
#' @param density.params list of parameters passed to [density.spherical]
#' @param nlevels integer. Number of contour levels for plotting
# #' @param type character. Type of plot: `'contour'` for contour lines,
# #' `'contour_filled'` for filled contours, or `'image'` for a raster image.
#' @param add logical. Whether the contours should be added to an existing plot.
#' @param col colour(s) for the contour lines drawn. If `NULL`, lines are color based on `col.palette`.
#' @param col.palette a color palette function to be used to assign colors in the plot.
#' @param col.params list. Arguments passed to `col.palette`
#' @param ... optional parameters passed to [graphics::image()] or [graphics::contour()].
#'
#' @name stereo_contour
#' @aliases contour image density-plot density_plot
#'
#' @family stereo-plot
#' @seealso [count_points()], [density]
#' @importFrom graphics image.default .filled.contour contour
#'
#' @returns list containing the stereographic x and coordinates of of the grid,
#' the counts, and the density.
#'
#' @examples
#' set.seed(20250411)
#' x <- rfb(100, mu = Line(120, 10), k = 5, A = diag(c(-1, 0, 1)))
#'
#' contour(x)
#'
#' x_density <- density(x)
#' contourf(x_density,
#'   col.params = list(direction = -1, begin = .05, end = .95, alpha = .75)
#' )
#' stereo_point(x, col = "black", pch = 21)
#'
#' image(x)
#' stereo_point(x, col = "lightgrey", pch = 21)
#'
#'
#' # complete example:
#' par(mfrow = c(1, 2))
#' wp <- 6 / ifelse(is.na(example_planes_df$quality), 6, example_planes_df$quality)
#' my_planes <- Plane(example_planes_df$dipdir, example_planes_df$dip)
#' fabric_p <- vollmer(my_planes)["D"]
#' my_planes_eig <- ot_eigen(my_planes)
#'
#' stereoplot(guides = TRUE, col = "grey96")
#' points(my_planes, col = "grey", pch = 16, cex = .5)
#' contour(my_planes, add = TRUE, weights = wp)
#' points(my_planes_eig$vectors[3, ], col = "black", pch = 16)
#' lines(my_planes_eig$vectors[3, ], ang = 90, col = "black", pch = 16)
#' title(
#'   main = "Planes",
#'   sub = paste0(
#'     "N: ", nrow(my_planes), " | Fabric strength: ", round(fabric_p, 2),
#'     "\nLambert equal area, lower hemisphere projection"
#'   )
#' )
#'
#' my_lines <- Line(example_lines_df$trend, example_lines_df$plunge)
#' wl <- 6 / ifelse(is.na(example_lines_df$quality), 6, example_lines_df$quality)
#' fabric_l <- vollmer(my_lines)["D"]
#'
#' stereoplot(guides = TRUE, col = "grey96")
#' points(my_lines, col = "grey", pch = 16, cex = .5)
#' contour(my_lines, add = TRUE, weights = wl)
#' points(sph_mean(my_lines, w = wl), col = "black", pch = 16)
#' title(
#'   main = "Lines",
#'   sub = paste0(
#'     "N: ", nrow(my_lines), " | Fabric strength: ", round(fabric_l, 2),
#'     "\nLambert equal area, lower hemisphere projection"
#'   )
#' )
NULL

#' @keywords internal
#' @rdname stereo_contour
#' @noRd
#' @importFrom viridis viridis
stereo_density <- function(x,
                           # kamb = TRUE, FUN = exponential_kamb, ngrid = 128L, sigma = 3,
                           # vmf_hw = NULL, vmf_optimal = c("cross", "rot"),
                           # weights = NULL, upper.hem = FALSE, r = 1,
                           density.params = list(),
                           type = c("contour", "contour_filled", "image"), nlevels = 10L,
                           col.palette = viridis::viridis, col = NULL, add = FALSE, col.params = list(),
                           ...) {
  type <- match.arg(type)

  if (inherits(x, "sph_density")) {
    d <- x
  } else {
    d <- do.call(density.spherical, append(list(x = x), density.params))
  }
  densities <- d$density

  if (isFALSE(add)) stereoplot(guides = FALSE)

  if (type == "image") {
    col.params <- append(list(n = nlevels), col.params)
    col <- do.call(col.palette, col.params)

    graphics::image.default(
      x = d$x, y = d$y,
      z = densities,
      col = col,
      asp = 1,
      axes = FALSE,
      frame.plot = FALSE,
      add = TRUE,
      ...
    )
  } else if (type == "contour") {
    if (is.null(col)) {
      levels <- pretty(range(densities, na.rm = TRUE), nlevels)
      col.params <- append(list(n = length(levels) - 1), col.params)
      col <- do.call(col.palette, col.params)
    }

    graphics::contour(
      x = d$x, y = d$y,
      z = densities,
      levels = levels,
      col = col,
      asp = 1,
      axes = FALSE,
      frame.plot = FALSE,
      add = TRUE,
      ...
    )
  } else {
    levels <- pretty(range(densities, na.rm = TRUE), nlevels)
    col.params <- append(list(n = length(levels) - 1), col.params)
    col <- do.call(col.palette, col.params)

    graphics::.filled.contour(
      x = d$x, y = d$y,
      z = densities,
      levels = levels,
      col = col
    )
  }

  invisible(densities)
}

#' @rdname stereo_contour
#' @exportS3Method stats::contour
contour.spherical <- function(x, add = FALSE, density.params = list(),
                              nlevels = 10L, col.palette = viridis::viridis, col = NULL,
                              col.params = list(), ...) {
  stereo_density(x,
    type = "contour",
    add = add,
    density.params = density.params,
    nlevels = nlevels,
    col = col,
    col.palette = col.palette,
    col.params = col.params,
    ...
  )
}

#' @rdname stereo_contour
#' @exportS3Method stats::contour
contour.sph_density <- function(x, add = FALSE, density.params = list(),
                                nlevels = 10L, col.palette = viridis::viridis, col = NULL,
                                col.params = list(), ...) {
  stereo_density(x,
    type = "contour",
    add = add,
    density.params = density.params,
    nlevels = nlevels,
    col = col,
    col.palette = col.palette,
    col.params = col.params,
    ...
  )
}


#' @rdname stereo_contour
#' @export
contourf <- function(x, add = FALSE, density.params = list(),
                     nlevels = 10L, col.palette = viridis::viridis,
                     col.params = list()) {
  stereo_density(x,
    type = "contour_filled",
    add = add,
    density.params = density.params,
    nlevels = nlevels,
    col.palette = col.palette,
    col.params = col.params
  )
}

#' @rdname stereo_contour
#' @exportS3Method stats::image
image.spherical <- function(x, add = FALSE, density.params = list(),
                            nlevels = 10L, col.palette = viridis::viridis,
                            col.params = list(), ...) {
  stereo_density(x,
    type = "image",
    add = add,
    density.params = density.params,
    nlevels = nlevels,
    col.palette = col.palette,
    col.params = col.params,
    ...
  )
}

#' @rdname stereo_contour
#' @exportS3Method stats::image
image.sph_density <- function(x, add = FALSE, density.params = list(),
                              nlevels = 10L, col.palette = viridis::viridis,
                              col.params = list(), ...) {
  stereo_density(x,
    type = "image",
    add = add,
    density.params = density.params,
    nlevels = nlevels,
    col.palette = col.palette,
    col.params = col.params,
    ...
  )
}
