#' Prepare points and lines for ggplot
#'
#' @inheritParams plot.spherical
#' @param ... [`<data-masking>`][rlang::args_data_masking()] Name-value pairs. The name gives the name of the column in the output.
#' The value can be:
#' \itemize{
#'  \item A vector of length 1, which will be recycled to the correct length.
#'  \item A vector the same length as the current group (or the whole data frame if ungrouped).
#'  \item `NULL`, to remove the column.
#'  \item A data frame or tibble, to create multiple columns in the output.
#' }
#' @param d numeric. Cone angle (small circle radius) in degrees. `90` (the default) produces great circles.
#' @param n integer. Resolution of line.
#'
#' @import dplyr
#'
#' @return data.frame
#'
#' @name prepare-ggplot
#'
#' @examples
#' if (require("mapproj")) {
#'   x <- Plane(120, 85)
#'   ggstereo() +
#'     ggplot2::geom_point(data = gg(x), ggplot2::aes(x, y), color = "red") +
#'     ggplot2::geom_path(data = ggl(x), ggplot2::aes(x, y), color = "red")
#'
#'   x2 <- Line(120, 5)
#'   ggstereo() +
#'     ggplot2::geom_point(data = gg(x2), ggplot2::aes(x, y), color = "darkgreen") +
#'     ggplot2::geom_path(
#'       data = ggl(x2, d = 8),
#'       ggplot2::aes(x, y, group = group), color = "darkgreen"
#'     )
#'
#'   x3 <- Plane(137, 71)
#'   ggstereo() +
#'     ggplot2::geom_point(data = gg(x3), ggplot2::aes(x, y), color = "darkgreen") +
#'     ggplot2::geom_path(
#'       data = ggl(x3, d = 90),
#'       ggplot2::aes(x, y, group = group), color = "darkgreen", lwd = 1
#'     ) +
#'     ggplot2::geom_path(
#'       data = ggl(x3, d = 90 + 11),
#'       ggplot2::aes(x, y, group = group, color = "sde <90")
#'     ) +
#'     ggplot2::geom_path(
#'       data = ggl(x3, d = 90 - 11),
#'       ggplot2::aes(x, y, group = group, color = "sde >90")
#'     )
#' }
NULL

#' @rdname prepare-ggplot
#' @export
gg <- function(x, ...) {
  stopifnot(is.spherical(x))
  azi <- inc <- NULL
  if (is.Plane(x) | is.Pair(x)) {
    x[, 1] <- (180 + x[, 1]) %% 360
    x[, 2] <- 90 - x[, 2]
  } else {
    x <- Line(x)
  }

  xdf <- data.frame(cbind(x))
  colnames(xdf) <- c("azi", "inc")

  xdf$x <- 180 - xdf$azi
  xdf$y <- xdf$inc

  dplyr::bind_cols(xdf, ...)
}


#' @rdname prepare-ggplot
#' @export
ggl <- function(x, ..., d = 90, n = 1e3) {
  id <- NULL
  stopifnot(is.spherical(x))

  if (n %% 2 > 0) n <- n + 1

  if (is.Plane(x) | is.Fault(x)) {
    # x[, 1] <- 180 + x[, 1]
    x[, 2] <- 90 - x[, 2]
  } else {
    x <- Line(x)
  }

  nx <- nrow(x)
  if (length(d) == 1) d <- rep(d, nx)

  xdf <- data.frame(cbind(x), ...) |>
    dplyr::mutate(id = dplyr::row_number(), d = d)

  zaxis <- Vec3(0, 0, 1)

  res <- matrix(ncol = 3, nrow = n * nx) |>
    as.data.frame()
  colnames(res) <- c("x", "y", "id")

  for (i in seq_len(nx)) {
    D <- Line(
      seq(0, 360, length.out = n),
      rep(90 - d[i], n)
    ) |>
      Vec3()

    strike <- as.Line(x[i, ])[1,1] - 90

    rotaxis <- Line(strike, 0) |> Vec3()
    D1 <- rotate(D, zaxis, deg2rad(strike))
    rotangle <- angle(zaxis, as.Line(x[i, ]))

    if (d[i] < 90 & is.Plane(x)) {
      d[i] <- 180 - d[i]
      D1 <- -D1
      x[i, 1] <- x[i, 1] + 180
      x[i, 2] <- 90 - x[i, 2]
    }

    # if (d[i] < 90) {
    #   k <- -1
    # } else {
    #   k <- 1
    # }
    k <- if (d[i] < 90) -1 else 1

    # lower hemisphere
    D_rot <- rotate(D1, rotaxis, k * rotangle) |> Line()
    D_fixed <- .fix_inc(az = D_rot[, 1], inc = D_rot[, 2])

    if (d[i] != 90 & d[i] > as.Line(x[i, ])[1,2]) {
      D_rotrot_fixed <- D_fixed

      dangle <- angle(Line(D_fixed[, 1], D_fixed[, 2]), as.Line(x[i, ]))
      cond <- dplyr::near(d[i], dangle)

      D_fixed[cond, 1] <- D_rotrot_fixed[cond, 1]
      D_fixed[cond, 2] <- D_rotrot_fixed[cond, 2]
    }

    D_fixed2 <- data.frame(x = (180 - D_fixed[, 1]), y = D_fixed[, 2], id = i)
    D_fixed3 <- if (d[i] < 90) D_fixed2 else utils::tail(D_fixed2, n = n / 2)

    idx <- ((i - 1) * n + 1):(i * n)
    res[idx, ] <- D_fixed3
  }

  res$group <- as.character(res$id)
  merge(res, xdf, by = "id")
}

#' Stereoplot Perimeter
#'
#' Adds a frame to the stereographic projection
#'
#' @inheritParams gg
#' @param color,fill,lwd Graphical parameters
#' @param ... optional graphical parameters passed to [ggplot2::geom_polygon()]
#'
#' @export
#' @importFrom ggplot2 aes geom_polygon
ggframe <- function(n = 1e4, color = "black", fill = NA, lwd = 1, ...) {
  prim.lat <- rep(c(0), times = n)
  prim.l1 <- seq(0, 180, length = n / 2)
  prim.l2 <- seq(-180, 0, length = n / 2)
  prim.long <- c(prim.l1, prim.l2)

  prim_df <- data.frame(prim.long, prim.lat)

  geom_polygon(aes(x = prim.long, y = prim.lat), data = prim_df, color = color, fill = fill, lwd = lwd, ..., inherit.aes = FALSE)
}

#' @importFrom ggplot2 aes geom_path
ggstereo_grid <- function(d = 10, rot = 0, ...) {
  x <- y <- group <- NULL
  # small circles
  sm1 <- seq(d, 90, d)
  sm2 <- seq(0, 90 - d, d)
  sm <- c(sm1, sm2)
  np <- Line(c(rep(0, length(sm1)), rep(180, length(sm2))), rep(0, length(sm)))

  # great circles
  dips <- seq(-90 + d, 90 - d, d)
  ep <- Plane((rep(90, length(dips)) * sign(dips)) %% 360, abs(dips))

  zp <- Plane(c(90, 0), c(90, 90))

  if (rot != 0) {
    zaxis <- Line(0, 90)
    np <- rotate(np, zaxis, rot)
    ep <- rotate(ep, zaxis, rot)
    zp <- rotate(zp, zaxis, rot)
  }

  sm_ggl <- ggl(np, d = sm)
  gc_ggl <- ggl(ep)
  zp_ggl <- ggl(zp)


  geom_path(data = dplyr::bind_rows(sm_ggl, gc_ggl, zp_ggl), mapping = aes(x, y, group = group), ..., inherit.aes = FALSE)
}


#' Stereonet using ggplot
#'
#' @param data Default dataset to use for plot. If not already a data.frame,
#' will be converted to one by [ggplot2::fortify()]. If not specified, must be
#' supplied in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot. If not
#' specified, must be supplied in each layer added to the plot.
#' @param earea logical. Whether the projection is equal-area ("Schmidt net")
#' (`TRUE`, the default), or equal-angle ("Wulff net") (`FALSE`).
#' @param grid.spacing numeric. Grid spacing in degree
#' @param grid.rot numeric. Angle (in degrees) to rotate the grid.
#' @param centercross logical. Whether a center cross should be added.
#' @param grid logical. Whether a gid should be added.
#' @param ... argument passed to [ggplot2::geom_polygon()]
#'
#' @import ggplot2
#' @importFrom ggplot2 aes annotate coord_map element_blank element_text ggplot scale_x_continuous scale_y_continuous theme
#'
#' @return ggplot
#' @export
#'
#' @examples
#' if (require("mapproj")) {
#'   test_data <- rbind(
#'     rvmf(100, mu = Line(90, 45), k = 10),
#'     rvmf(50, mu = Line(0, 0), k = 20)
#'   ) |> as.line()
#'
#'   ggstereo(grid = TRUE) +
#'     ggplot2::geom_point(data = gg(test_data), ggplot2::aes(x = x, y = y))
#'
#'   ggstereo(earea = FALSE, centercross = TRUE) +
#'     ggplot2::geom_point(data = gg(test_data), ggplot2::aes(x = x, y = y))
#' }
ggstereo <- function(data = NULL, mapping = aes(), earea = TRUE, centercross = TRUE, grid = FALSE, grid.spacing = 10, grid.rot = 0, ...) {
  # if(earea){
  #   crs = "+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0"
  # } else {
  #   crs = "+proj=stere +lat_0=90 +lon_0=0 +x_0=0 +y_0=0"
  # }
  rlang::check_installed("mapproj", reason = "to use `coord_map()`")

  ggplot(data = data, mapping = mapping) +
    # theme_void() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      # legend.title = element_blank()
    ) +
    {
      if (grid) {
        ggstereo_grid(d = grid.spacing, rot = grid.rot, color = "grey90", lwd = .2)
      }
    } +
    ggframe(...) +
    annotate("point", x = 0, y = 90, pch = as.numeric(centercross) * 3) +
    scale_y_continuous(limits = c(0, 90)) +
    scale_x_continuous(limits = c(-180, 180)) +
    coord_map(ifelse(earea, "azequalarea", "stereographic"), orientation = c(90, 0, 0))
  # coord_sf(crs = crs, default_crs = crs)
}

ignore_unused_imports <- function() {
  mapproj::mapproject
}




.full_hem <- function(azi, inc) {
  inc <- inc + 90

  # Create upper hemisphere reflection
  azi2 <- (azi + 180) %% 360
  inc2 <- 180 - inc

  # Combine and return as matrix
  azi3 <- c(azi, azi2)
  inc3 <- c(inc, inc2)

  cbind("inc" = inc3, "azi" = azi3)
}


#' Stereonet contouring using ggplot
#'
#' @inheritParams ggstereo
#' @param ngrid integer. Resolution of density calculation.
#' @param hw numeric. Kernel bandwidth in degree.
#' @param optimal_bw character. Calculates an optimal kernel bandwidth
#' using the cross-validation algorithm (`'cross'`) or the rule-of-thumb (`'rot'`)
#' suggested by Garcia-Portugues (2013). Ignored when `hw` is specified.
#' @param norm logical. Should the densities be normalized?
#' @param smooth logical. Whether [ggplot2::geom_tile()] should be used for plotting.
#' @param threshold numeric. Cut-off for low densities.
#' @param ... arguments passed to [ggplot2::geom_contour()], [ggplot2::geom_contour_filled()], or [ggplot2::geom_tile()]
#'
#' @return ggplot
#'
#' @references Garcia Portugues, E. (2013). Exact risk improvement of
#' bandwidth selectors for kernel density estimation with directional data.
#' Electronic Journal of Statistics, 7, 1655-1685.
#'
#' @import ggplot2
#' @importFrom Directional vmf.kerncontour euclid vmfkde.tune
#' @name ggstereocontour
#' @examples
#' if (require("mapproj")) {
#'   test_data <- rbind.spherical(
#'     rvmf(100, mu = Line(90, 45), k = 10),
#'     rvmf(50, mu = Line(0, 0), k = 20)
#'   )
#'
#'   ggstereo() +
#'     geom_contourf_stereo(gg(test_data)) +
#'     ggplot2::scale_fill_viridis_d(option = "A") +
#'     # guides(fill = guide_colorsteps(barheight = unit(8, "cm"), show.limits = TRUE)) +
#'     geom_contour_stereo(gg(test_data), color = "grey") +
#'     ggplot2::geom_point(data = gg(test_data), ggplot2::aes(x = x, y = y), color = "lightgrey") +
#'     ggframe()
#'
#'   ggstereo() +
#'     geom_contourf_stereo(gg(test_data), norm = TRUE, bins = 50, threshold = .1) +
#'     ggplot2::scale_fill_viridis_d(option = "A")
#' }
NULL

#' @rdname ggstereocontour
#' @export
#' @importFrom ggplot2 aes geom_contour
geom_contour_stereo <- function(data, ngrid = 200, hw = NULL, optimal_bw = c("cross", "rot"), norm = FALSE, threshold = 0, ...) {
  Long <- Lat <- Density <- NULL
  xtot <- .full_hem(azi = data$azi, inc = data$inc)

  dens <- vmf_kerncontour(xtot, hw = hw, kernel_method = optimal_bw, ngrid = ngrid)
  res <- expand.grid(Lat = dens$lat - 90, Long = dens$long - 180)
  res$Density <- c(dens$den)
  if (norm) {
    res$Density <- normalize(res$Density)
  }
  res$Density[res$Density <= threshold] <- NA

  geom_contour(data = res, aes(x = -Long, y = Lat, z = Density), ...)
}


#' @rdname ggstereocontour
#' @export
#' @importFrom ggplot2 aes geom_contour_filled geom_tile
geom_contourf_stereo <- function(data, ngrid = 200, hw = NULL, optimal_bw = c("cross", "rot"), norm = FALSE, smooth = FALSE, threshold = 0, ...) {
  Long <- Lat <- Density <- NULL
  xtot <- .full_hem(azi = data$azi, inc = data$inc)

  dens <- vmf_kerncontour(xtot, hw = hw, kernel_method = optimal_bw, ngrid = ifelse(smooth, 3 * ngrid, ngrid))
  res <- expand.grid(Lat = dens$lat - 90, Long = dens$long - 180)
  res$Density <- c(dens$den)
  if (norm) {
    res$Density <- normalize(res$Density)
  }
  res$Density[res$Density <= threshold] <- NA

  if (smooth) {
    geom_tile(data = res, aes(x = -Long, y = Lat, fill = Density), ...)
  } else {
    geom_contour_filled(data = res, aes(x = -Long, y = Lat, z = Density), ...)
  }
}

#
# opening.angle <- function(x, oa.max, oa.min, oa.bw = 10) {
#   sub1 <- subset(x, inc <= oa.max) # max degree increment to include
#   sub2 <- subset(sub1, inc >= oa.min) # min degree increment to include
#   peakx <- density(sub2$az, bw = oa.bw)$x[which(diff(sign(diff(density(sub2$az, bw = oa.bw)$y))) == -2)] # determine the x value of all data peaks
#   peaky <- density(sub2$az, bw = oa.bw)$y[which(diff(sign(diff(density(sub2$az, bw = oa.bw)$y))) == -2)] # determine the y value of all data peaks
#   if (length(peaky) > 4) {
#     peakn <- length(peaky)
#     for (peakn in 1:peakn - 4) {
#       peakx <- peakx[-(which.min(peaky))]
#       peaky <- peaky[-(which.min(peaky))]
#       if (length(peaky) == 4) break
#     }
#     OA <- round((((peakx[3] - peakx[2]) + (360 - peakx[4] + peakx[1])) / 2), digits = 0)
#   } else {
#     OA <- round((((peakx[3] - peakx[2]) + (360 - peakx[4] + peakx[1])) / 2), digits = 0)
#   }
#   def.temp <- round(OA.temp(OA), digits = 0)
#   out <- list(peakx, peaky, OA, sub2$az, def.temp)
# }
