#' Spatial interpolation
#' 
#' Inverse distance weighted spatial interpolation of plane or line objects.
#'
#' @param x numeric vector, array, or object of class `"line"` or `"plane"`
#' @param coords a `"sf"` object containing the geographic coordinates of `x` measurements
#' @param grid (optional) Point object of class `"sf"`.`
#' @param lon_range,lat_range (optional) Numeric vector specifying the minimum and maximum longitudes and latitudes (ignored if `grid` is specified).
#' @param gridsize 	Numeric. Target spacing of the regular grid in decimal degree. Default is `2.5`. (is ignored if grid is specified)
#' @param min_data Integer. Minimum number of data per kernel. Default is `3`
#' @param threshold Numeric. Threshold for deviation of direction. Default is `25`
#' @param arte_thres Numeric. Maximum distance (in km) of the grid point to the next data point. Default is `200`
#' @param dist_weight Distance weighting method which should be used: `"linear"`, or `"inverse"` (the default).
#' @param idp Numeric. The weighting power of inverse distance. When set to `0`, no weighting is applied.
#' @param dist_threshold Numeric. Distance weight to prevent overweight of data nearby (0 to 1). Default is `0.1`
#' @param R Numeric value or vector specifying the kernel half-width, i.e. the search radius (in km). Default is `1`
#' 
#' @returns list
#' @importFrom tectonicr dist_greatcircle  
#' @importFrom sf st_transform st_coordinates st_is st_as_sf st_make_grid st_crs st_bbox
#' @details Based on [tectonicr::stress2grid()]
#' 
#' @seealso [tectonicr::stress2grid()], [v_mean()], [v_sd()]
#' @export
#' @examples
#' \dontrun{
#' data <- read_strabo_JSON("E:/Lakehead/Field work/StraboSpot_07_02_2023.json", dataset = "TS")
#' ps <- data$data |> dplyr::mutate(dipdir = (strike + 90) %% 360) |> dplyr::filter(type == "planar_orientation" & !(feature_type %in% c("other", "vector", "option_13")))
#' ps_vec <- structr::as.plane(cbind(ps$dipdir, ps$dip))
#' structr::spatial_interpolation(x=ps_vec, coords = ps, gridsize = .05, R = 15, dist_threshold = 0.01, threshold = Inf)
#' }
spatial_interpolation <- function(x,
                                  coords,
                                  grid = NULL,
                                  lon_range = NULL,
                                  lat_range = NULL,
                                  gridsize = .1,
                                  min_data = 3,
                                  threshold = Inf,
                                  arte_thres = 100,
                                  dist_weight = c("inverse", "linear"),
                                  idp = 1.0,
                                  dist_threshold = 0.01,
                                  R = 10,
                                  ...) {
  # stopifnot(inherits(coords, "sf"), is.numeric(threshold), is.numeric(arte_thres),
  #   arte_thres > 0, is.numeric(dist_threshold), is.numeric(R), is.numeric(idp), 
  # )

  transform <- FALSE
  if (is.spherical(x)) {
    v <- to_vec(x)
    transform <- TRUE
  } else {
    v <- vec2mat(x)
  }

  min_data <- as.integer(ceiling(min_data))
  dist_weight <- match.arg(dist_weight)

  # pre-allocating
  #N <- nrow(x)
  lat <- lon <- numeric()

  coords <- coords |>
    sf::st_transform(crs = "WGS84") |>
    sf::st_coordinates()

  datas <- cbind(
    lon = coords[, 1],
    lat = coords[, 2],
    x = v[, 1],
    y = v[, 2],
    z = v[, 3]
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
  G <- sf::st_coordinates(grid)

  N <- numeric(nrow(G))

  SH <- c()
  for (i in seq_along(G[, 1])) {
    distij <- tectonicr::dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1])

    if (min(distij) <= arte_thres) {
      # for (k in seq_along(R)) {
      # R_search <- R
      ids_R <-
        which(distij <= R) # select those that are in search radius

      N_in_R <- length(ids_R)

      if (N_in_R < min_data) {
        # not enough data within search radius
        sd_vec <- NA
        mean_vec <- cbind(NA, NA, NA)
        mdr <- NA
      } else if (N_in_R == 1) {
        sd_vec <- 0
        mean_vec <- datas[ids_R, 3:5]
        mdr <- distij[ids_R] / R
      } else {
        mdr <- mean(distij[ids_R], na.rm = TRUE) / R
        dist_threshold_scal <- R * dist_threshold

        if (dist_weight == "linear") {
          w <- R + 1 - max(dist_threshold_scal, distij[ids_R])
        } else {
          w <- 1 / (max(dist_threshold_scal, distij[ids_R]))^idp
        }
        
        # mean value
        mean_vec <- v_mean(datas[ids_R, 3:5], w)
        sd_vec <- v_delta(datas[ids_R, 3:5], w)
      }
      SH.ik <- c(
        lon = G[i, 1],
        lat = G[i, 2],
        x = mean_vec[1],
        y = mean_vec[2],
        z = mean_vec[3],
        sd = tectonicr::rad2deg(sd_vec),
        mdr = mdr,
        N = N_in_R
      )

      if (SH.ik['sd'] <= threshold & !is.na(SH.ik['sd'])) {
        SH <- rbind(SH, SH.ik)
      }
      # }
    }
  }

  res_coords <- data.frame(X = SH[, 1], Y = SH[, 2])

  if (transform) {
    vec <- to_spherical(cbind(x = SH[, 3], y = SH[, 4], z = SH[, 5]), class(x))
  } else {
    vec <- cbind(x = SH[, 3], y = SH[, 4], z = SH[, 5])
  }

  stats <- data.frame(
    sd = SH[, "sd"],
    mdr = SH[, "mdr"],
    N = SH[, "N"]
  )

  res <- list(
    coords = res_coords,
    mean = vec,
    stats = stats
  )

  return(res)
}
