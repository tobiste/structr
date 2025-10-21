#' Spatial interpolation
#'
#' Inverse distance weighted spatial interpolation of plane or line objects.
#'
#' @param x numeric vector, array, or object of class `"line"` or `"plane"`
#' @param coords a `"sf"` object containing the geographic coordinates of `x`
#' measurements
#' @param grid (optional) Point object of class `"sf"`.
#' @param gridsize 	Numeric. Target spacing of the regular grid in decimal
#'  degree. Default is `2.5`. Iignored if grid is specified.
#' @param min_data Integer. Minimum number of data per kernel. Default is `3`
#' @param threshold Numeric. Threshold for deviation of direction. Default is `25`
#' @param arte_thres Numeric. Maximum distance (in km) of the grid point to the
#'  next data point. Default is `200`
#' @param dist_weighting Distance weighting method which should be used:
#'  `"linear"`, or `"inverse"` (the default).
#' @param idp Numeric. The weighting power of inverse distance. When set to `0`,
#'  no weighting is applied.
#' @param dist_threshold Numeric. Distance weight to prevent overweight of data
#'  nearby (0 to 1). Default is `0.1`
#' @param R_range Numeric value or vector specifying the kernel half-width, i.e.
#'  the search radius (in km). Default is `1`
#' @param .compact logical. Run [compact_grid()] at the end
#' @param lon_range,lat_range two column vector specifying the coordinate range.
#' Ignored when grid is specified.
#'
#' @returns list
#' @importFrom tectonicr dist_greatcircle
#' @importFrom dplyr bind_rows mutate select
#' @importFrom sf st_transform st_coordinates st_is st_as_sf st_make_grid st_crs st_bbox
#' @details Based on [tectonicr::stress2grid()]
#'
#' @seealso [tectonicr::stress2grid()], [sph_mean()], [delta()]
#' @noRd
#' @examples
#' ps_vec <- rvmf() |> Line()
#' ps <- data.frame(x = stats::runif(100, 50, 60), y = stats::runif(100, 40, 45)) |>
#'   st_as_sf(coords = c("x", "y"))
#'
#' spatial_interpolation(x = ps_vec, coords = ps, gridsize = .5, .compact = FALSE)
spatial_interpolation <- function(x,
                                  coords,
                                  grid = NULL,
                                  lon_range = NULL,
                                  lat_range = NULL,
                                  gridsize = 1L,
                                  min_data = 3L,
                                  max_data = Inf,
                                  max_sd = Inf,
                                  min_dist_threshold = Inf,
                                  dist_weighting = c("inverse", "linear", "none"),
                                  idp = 1,
                                  dist_threshold = 0.1,
                                  R_range = seq(1, 10, 1),
                                  .compact = TRUE) {
  stopifnot(
    is.spherical(x),
    is.numeric(gridsize) && length(gridsize) == 1,
    is.numeric(max_sd) | is.infinite(max_sd),
    is.numeric(max_data) | is.infinite(max_data),
    is.numeric(min_data) | is.infinite(min_data),
    max_data >= min_data,
    is.numeric(min_dist_threshold),
    is.numeric(dist_threshold),
    min_dist_threshold > 0 && length(min_dist_threshold) == 1,
    is.numeric(R_range),
    is.numeric(idp) && length(idp) == 1
  )

  v <- Line(x) |> unclass()

  arte_thres <- min_dist_threshold
  threshold <- max_sd
  min_data <- as.integer(ceiling(min_data))
  dist_weighting <- match.arg(dist_weighting)
  w_distance_fun <- if (dist_weighting == "linear") tectonicr:::dist_weighting_linear else tectonicr:::dist_weighting_inverse

  # colnames_x <- colnames(x)

  # pre-allocating
  n_vecs <- nrow(x)
  lat <- lon <- numeric(n_vecs)
  N <- md <- R <- numeric()

  if (dist_weighting == "none") idp <- 0

  x_coords <- sf::st_coordinates(coords)
  datas <- cbind(
    lon = x_coords[, 1],
    lat = x_coords[, 2],
    azimuth = v[, 1],
    plunge = v[, 2]
  ) #|> as.matrix()

  if (is.null(grid)) {
    # Regular grid
    if (is.null(lon_range) | is.null(lat_range)) {
      lon_range <- range(datas[, 1], na.rm = TRUE)
      lat_range <- range(datas[, 2], na.rm = TRUE)
    }

    grid <- sf::st_bbox(
      c(
        xmin = lon_range[1],
        xmax = lon_range[2],
        ymin = lat_range[1],
        ymax = lat_range[2]
      ),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(
        cellsize = gridsize,
        what = "centers",
        offset = c(lon_range[1], lat_range[1])
      ) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- unname(sf::st_coordinates(grid))
  R_seq <- seq_along(R_range)

  res <- lapply(seq_along(G[, 1]), function(i) {
    # for(i in seq_along(G[, 1])){
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1])
    if (max_data < Inf) distij <- distij[tectonicr:::which.nsmallest(distij, max_data)] # select the `max_data` nearest locations

    if (min(distij) <= min_dist_threshold) {
      t(vapply(R_seq, function(k) {
        # for(k in R_seq){
        R_search <- R_range[k]
        ids_R <- (distij <= R_search) # select those that are in search radius
        N_in_R <- sum(ids_R)

        if (N_in_R < min_data) {
          # not enough data within search radius
          sd_vec <- NA
          mean_vec <- cbind(NA, NA)
          md <- NA
        } else if (N_in_R == 1) {
          sd_vec <- 0
          mean_vec <- datas[ids_R, 3:4]
          md <- distij[ids_R]
        } else {
          md <- mean(distij[ids_R], na.rm = TRUE)

          # distance weighting
          w <- w_distance_fun(R_search, dist_threshold, distij[ids_R], idp)


          # mean vector and spherical standard deviation
          vecs <- datas[ids_R, 3:4] |> as.Line()
          mean_vec <- sph_mean(vecs, w) |> unclass()
          sd_vec <- delta(vecs, w)
        }
        c(
          lon = G[i, 1],
          lat = G[i, 2],
          azimuth = mean_vec[1, 1],
          plunge = mean_vec[1, 2],
          # z = mean_vec[1,3],
          delta = sd_vec,
          R = R_search,
          md = md,
          N = N_in_R
        )
      }, FUN.VALUE = numeric(8)))
    }
  }) |>
    lapply(as.data.frame) |>
    dplyr::bind_rows() |>
    dplyr::mutate(mdr = md / R, N = as.integer(N)) |>
    dplyr::select(-md) |>
    # dplyr::filter(delta <= threshold, !is.na(delta)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)

  # vec <- Line(res$azimuth, res$plunge)
  #
  # if (is.Plane(x)) {
  #   res <- res |> select(-c("azimuth", "plunge"))
  #   vec <- Plane(vec)
  #   res$dipdir <- vec[, 1]
  #   res$dip <- vec[, 2]
  # } else if (is.Vec3(x)) {
  #   res <- res |> select(-c("azimuth", "plunge"))
  #   res$x <- vec[, 1]
  #   res$y <- vec[, 2]
  #   res$z <- vec[, 2]
  # }

  if (.compact) res <- compact_grid(res)

  return(res)
}

#' Compact grid
#'
#' Filter spatial interpolation containing a range of search radii or kernel
#' half widths to find smallest wavelength (R) with the least spherical sd.
#'
#' @param grid output of [tectonicr::stress2grid()],
#' [tectonicr::PoR_stress2grid()], or [tectonicr::kernel_dispersion()]
#' @returns \code{sf} object
#'
#' @importFrom stats aggregate
#'
#' @noRd
compact_grid <- function(grid) {
  lon <- lat <- x <- y <- z <- dipdir <- dip <- R <- numeric()
  group <- character()

  data <- subset(grid, !is.na(grid$x))
  data$group <- paste(lon, lat)

  temp <- aggregate(R ~ group, data, min, na.rm = TRUE)
  temp2 <- merge(temp, data, by.x = "group", by.y = "R")
  temp2$group <- NULL

  sf::st_as_sf(temp2)
}
