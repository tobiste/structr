#' Import orientation data from *StraboSpot*
#'
#' Reads the `XLS` format export of field book data, the `JSON` project file, or
#' the `txt` export of StraboMobile data from \url{strabospot.org/my_data} and
#' creates a list with the metadata, and the line or plane orientations.
#'
#' @param file the name of the file which the data are to be read from.
#' @param dataset character. name of the dataset extracted from `JSON` project
#' file.
#' @param tag_cols logical. Whether the Tag columns should be summarized in a
#' single column (may lead to duplicate rows).
#' @param sf logical. Whether the output should be a spatial `"sf"` object
#' using the Longitude and Latitude columns.
#' @importFrom readxl read_xlsx
#' @importFrom rjson fromJSON
#' @importFrom lubridate as_datetime
#' @importFrom dplyr rename mutate filter select full_join arrange pull bind_rows
#' @importFrom plyr rbind.fill
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect starts_with
#' @importFrom tibble tibble as_tibble
#' @importFrom sf st_as_sf
#' @returns `list` containing the following objects:
#' \describe{
#' \item{`data`}{`"tbl_df"` object or `"sf"` if `sf == TRUE`. Metadata.}
#' \item{`tags`}{`"tbl_df"` object. Tags and their descriptsions.}
#' \item{`planar`}{Plane elements. Same row IDs as in `data`.}
#' \item{`linear`}{Line elements. Same row IDs as in `data`.}
#' }
#' @name strabo
#' @examples
#' \dontrun{
#' file <- "E:/Lakehead/Field Work/StraboSpot_Output_11_01_2022_TS.xlsx"
#' dt <- read_strabo_xls(file)
#' stereoplot()
#' stereo_greatcircle(dt$planar, col = "lightgrey")
#' stereo_point(dt$linear)
#'
#' read_strabo_mobile("C:/Users/tobis/Downloads/StraboSpot_Search_06_16_2023.txt")
#'
#' read_strabo_JSON("G:/My Drive/Moss_Lake/data.json")
#' }
NULL

#' @rdname strabo
#' @export
read_strabo_xls <- function(file, tag_cols = FALSE, sf = TRUE) {
  Date <- Planar.Orientation.Dipdirection <- Planar.Orientation.Strike <-
    Linear.Sense <- Planar.Orientation.Movement <- temp <- Tag <- Linear.Orientation.Trend <-
    Planar.Orientation.Fault.Or.Sz.Type <- Planar.Orientation.Unix.Timestamp <- Linear.Orientation.Unix.Timestamp <-
    Planar.Orientation.Dip <- Linear.Orientation.Plunge <- Longitude <- Latitude <- NULL

  data <- readxl::read_xlsx(file, sheet = 1, skip = 2)
  colnames(data) <- make.names(colnames(data))

  data <- data |>
    dplyr::mutate(
      Date = lubridate::as_datetime(Date),
      Planar.Orientation.Dipdirection = (Planar.Orientation.Strike + 90) %% 360,
      Linear.Sense = ifelse(Planar.Orientation.Movement %in% c("left_lateral"), -1, NA),
      Linear.Sense = ifelse(Planar.Orientation.Movement %in% c("right_lateral"), 1, Linear.Sense),
      Linear.Sense = ifelse(Planar.Orientation.Fault.Or.Sz.Type %in% c("sinistral", "reverse"), -1, Linear.Sense),
      Linear.Sense = ifelse(Planar.Orientation.Fault.Or.Sz.Type %in% c("dextral", "normal", "dextral_normal"), 1, Linear.Sense)
    )
  if (tag_cols) {
    data <- data |>
      tidyr::pivot_longer(cols = tidyselect::starts_with("Tag"), names_to = "Tag", values_to = "temp") |>
      dplyr::filter(temp == "X") |>
      dplyr::select(-temp) |>
      dplyr::mutate(Tag = gsub("Tag:", "", Tag))
  }

  data.lines0 <- data |>
    dplyr::select(tidyselect::starts_with("Linear")) |>
    dplyr::filter(!is.na(Linear.Orientation.Trend))

  data.planes0 <- data |>
    dplyr::select(tidyselect::starts_with("Planar"), ) |>
    dplyr::filter(!is.na(Planar.Orientation.Dipdirection))


  data0 <- data |> dplyr::select(
    !tidyselect::starts_with(c("Planar", "Linear"))
  )
  data0$Planar.Orientation.Unix.Timestamp <- data$Planar.Orientation.Unix.Timestamp

  res <- dplyr::full_join(data.planes0, data.lines0,
    by = c("Planar.Orientation.Unix.Timestamp" = "Linear.Orientation.Unix.Timestamp"),
    multiple = "all"
  ) |>
    dplyr::left_join(data0,
      by = "Planar.Orientation.Unix.Timestamp",
      multiple = "all"
    ) |>
    dplyr::mutate(Linear.Orientation.Unix.Timestamp = NA) |>
    dplyr::select(colnames(data)) |>
    dplyr::select(-Linear.Orientation.Unix.Timestamp)

  planes <- as.plane(cbind(res$Planar.Orientation.Dipdirection, res$Planar.Orientation.Dip))
  lines <- as.line(cbind(res$Linear.Orientation.Trend, res$Linear.Orientation.Plunge))

  if (sf) {
    res <- sf::st_as_sf(res, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
  } else {
    res <- res
  }
  list(
    data = res,
    planar = planes,
    linear = lines
  )
}

#' @rdname strabo
#' @export
read_strabo_mobile <- function(file, sf = TRUE) {
  x <- Type <- No. <- Trd.Strk <- Plg.Dip <- Dipdir <- NULL
  data0 <- utils::read.table(x, header = TRUE, sep = "\t", colClasses = c("integer", rep("character", 3), rep("numeric", 11), "character"))

  lines0 <- dplyr::filter(data0, Type == "L")
  lines <- as.line(cbind(lines0$Trd.Strk, lines0$Plg.Dip))
  lines.meta <- lines0 |> dplyr::select(-c(No., Trd.Strk, Plg.Dip))
  rownames(lines) <- rownames(lines.meta) <- lines0$No.

  planes0 <- dplyr::filter(data0, Type == "P") |>
    mutate(Dipdir = (Trd.Strk + 90) %% 360)
  planes <- as.plane(cbind(planes0$Dipdir, planes0$Plg.Dip))
  planes.meta <- planes0 |> dplyr::select(-c(No., Dipdir, Trd.Strk, Plg.Dip))
  rownames(planes) <- rownames(planes.meta) <- planes0$No.


  if (sf) {
    lines.meta <- sf::st_as_sf(lines.meta, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
    planes.meta <- sf::st_as_sf(planes.meta, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
  }

  list(
    linear = lines,
    linear_data = lines.meta,
    planar = planes,
    planar_data = planes.meta
  )
}

#' @rdname strabo
#' @export
read_strabo_JSON <- function(file, sf = TRUE) {
  time <- tag_id <- tag_name <- NULL

  dat <- rjson::fromJSON(file = file)

  # read tags
  tags_list <- dat$project$project$tags
  tags_df <- spot_tags <- tibble::tibble()
  for (t in seq_along(tags_list)) {
    spots_t <- data.frame(spot = tags_list[[t]]$spots)
    tags_list[[t]]$spots <- NULL
    tags_list[[t]]$eon <- NULL

    df_t <- tibble::as_tibble(tags_list[[t]])
    tags_df <- dplyr::bind_rows(tags_df, df_t)
    if (nrow(spots_t) > 0) {
      spots_t$tag <- df_t$id
      spot_tags <- rbind(spot_tags, spots_t)
    }
  }
  colnames(tags_df) <- paste("tag", colnames(tags_df), sep = "_")

  spot_tags_df <- dplyr::left_join(
    unique(spot_tags),
    tags_df |> dplyr::select(tag_id, tag_name),
    by = c("tag" = "tag_id")
  ) |>
    dplyr::mutate(tag_name = paste0("tag:", tag_name)) |>
    # dplyr::select(-tag) |>
    tidyr::pivot_wider(names_from = "tag_name", values_from = "tag") |>
    dplyr::mutate(dplyr::across(!spot, function(x) ifelse(is.na(x), FALSE, TRUE))) |>
    unique() |>
    mutate(spot = as.character(spot)) |>
    rename(spot_id = spot)

  # datasets
  ds <- seq_along(dat$project$datasets)
  ds_name <- ds_name_id <- character()

  for (a in ds) {
    ds_name[a] <- dat$project$datasets[[a]]$name
    ds_name_id[a] <- names(dat$project$datasets)[[a]]
  }
  ds_df <- data.frame(ds_name = ds_name, ds_id = ds_name_id)

  # spot lists
  spots <- character()
  for (i in ds_name_id) {
    ds_spot_ids_i <- dat$project$datasets[[i]]$spotIds
    spots <- rbind(spots, cbind(ds_id = rep(i, length(ds_spot_ids_i)), spot_id = ds_spot_ids_i))
  }
  spots_df <- as.data.frame(spots) |> unique()

  # extract the specified dataset
  dsn <- data.frame(ds, ds_name) |>
    # dplyr::filter(ds_name == dataset) |>
    dplyr::pull(ds)
  ds_name_id <- names(dat$project$datasets)[dsn]


  # ds_spot_ids <- dat$project$datasets[[dsn]]$spotIds # spots in dataset
  # DSN <- dat$spotsDb[ds_spot_ids]
  # for (i in seq_along(dat$spotsDb)) {
  #   dat$spotsDb[[i]]
  # }


  spot_id <- spot <- notes <- date <- character()
  altitude <- longitude <- latitude <- gps_accuracy <- numeric()
  orient_df <- data.frame(check.names = FALSE)

  id <- numeric()
  for (i in seq_along(dat$spotsDb)) {
    spot_i <- dat$spotsDb[[i]]
    if (!is.null(spot_i$geometry)) {
      if (spot_i$geometry$type == "Point" & spot_i$type == "Feature") {
        spot_id[i] <- spot_i$properties$id
        spot[i] <- spot_i$properties$name
        date[i] <- spot_i$properties$date
        id[i] <- spot_i$properties$id

        notes_i <- spot_i$properties$notes
        notes[i] <- ifelse(is.null(notes_i), NA, notes_i)

        altitude_i <- spot_i$properties$altitude
        altitude[i] <- ifelse(is.null(altitude_i), NA, altitude_i)

        longitude[i] <- spot_i$geometry$coordinates[1]
        latitude[i] <- spot_i$geometry$coordinates[2]

        gps_accuracy_i <- spot_i$properties$gps_accuracy
        gps_accuracy[i] <- ifelse(is.null(gps_accuracy_i), NA, gps_accuracy_i)

        # orientation data
        orient_ls_i <- spot_i$properties$orientation_data
        if (!is.null(orient_ls_i)) {
          orient_df_i <- tibble::tibble()
          for (j in seq_along(orient_ls_i)) { # loop through each measurement
            orient_ls_ij <- orient_ls_i[[j]]

            # concatenate directional_indicators into one column
            if (length(orient_ls_ij$directional_indicators) > 1) {
              orient_ls_ij$directional_indicators <- toString(orient_ls_ij$directional_indicators)
            }

            # single P or L measurements
            if (is.null(orient_ls_ij$associated_orientation)) {
              # append all measurements of spot
              # orient_df_i <- plyr::rbind.fill(orient_df_i, as.data.frame(orient_ls_ij, check.names = FALSE)) |>
              orient_df_i <- plyr::rbind.fill(orient_df_i, as.data.frame(orient_ls_ij, check.names = FALSE))
              orient_df_i$spot <- spot[i]
              orient_df_i$spot_id <- spot_id[i]
              orient_df_i$associated <- FALSE
            } else {
              # P + L measurements

              L_ls <- orient_ls_ij$associated_orientation
              orient_ls_ij$associated_orientation <- NULL

              ## P measurement
              P <- orient_ls_ij |>
                as.data.frame(check.names = FALSE) |>
                unique()

              ## L measurement
              PL <- data.frame(check.names = FALSE)
              for (k in seq_along(L_ls)) {
                if (L_ls[[k]]$type != "planar_orientation") {
                  L <- L_ls[[k]] |>
                    as.data.frame(check.names = FALSE) |>
                    unique()
                  plunge <- L$plunge
                  trend <- L$trend
                  if (!is.na(plunge) & !is.na(trend)) {
                    L$unix_timestamp <- L$plunge <- L$trend <- NULL # don't need this column twice
                    colnames(L) <- paste("associated", colnames(L), sep = "_")
                    L$plunge <- plunge
                    L$trend <- trend

                    ## combine P + L measurements
                    PL_k <- cbind(P, L)
                    PL <- plyr::rbind.fill(PL, PL_k) |>
                      unique()
                  }
                  PL$associated <- TRUE
                }
              }

              # append all P+L measurements of spot
              if (nrow(orient_df_i) > 0) {
                orient_df_i <- plyr::rbind.fill(orient_df_i, PL) |>
                  unique()
                orient_df_i$spot <- spot[i]
                orient_df_i$spot_id <- spot_id[i]
              }

              # orient_ls_ij <- NULL
            }
            # combine all spots
            orient_df <- plyr::rbind.fill(orient_df, orient_df_i)
          }
        }
      }
    }
  }
  fieldbook <- tibble::tibble(spot_id = spot_id, spot = spot, longitude, latitude, altitude, note = notes, time = date, gps_accuracy, id) |>
    dplyr::mutate(time = lubridate::as_datetime(time)) |>
    dplyr::arrange(time) |>
    dplyr::distinct() |>
    dplyr::left_join(spots_df) |>
    dplyr::left_join(ds_df) |>
    dplyr::left_join(spot_tags_df) |>
    dplyr::select(-c(id))
  if (sf) fieldbook <- sf::st_as_sf(fieldbook, coords = c("longitude", "latitude"), remove = FALSE, crs = "WGS84", na.fail = FALSE)

  meta <- orient_df |>
    tibble::as_tibble() |>
    dplyr::distinct() |>
    dplyr::left_join(fieldbook, dplyr::join_by("spot_id")) |>
    dplyr::arrange(time) |>
    dplyr::select(-c(spot.y, ds_id, spot_id, unix_timestamp, modified_timestamp)) |>
    dplyr::rename(spot = spot.x)
  #
  if (sf) {
    meta <- sf::st_as_sf(meta)
  }

  planes <- as.plane(cbind(orient_df$strike + 90, orient_df$dip))
  lines <- as.line(cbind(orient_df$trend, orient_df$plunge))

  list(
    data = meta,
    tags = tags_df,
    planar = planes,
    linear = lines
  )
}
