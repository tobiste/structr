`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Import orientation data from *StraboSpot*
#'
#' Reads the `XLS` format export of field book data, the `JSON` project file, or
#' the `txt` export of StraboMobile data from \url{strabospot.org/my_data} and
#' creates a list with the metadata, and the line or plane orientations.
#'
#' @param file the name of the file which the data are to be read from.
#' @param tag_cols logical. Whether the Tag columns should be summarized in a
#' single column (may lead to duplicate rows).
#' @param sf logical. Whether the output should be a spatial `"sf"` object
#' using the Longitude and Latitude columns.
#' @importFrom readxl read_xlsx
#' @importFrom rjson fromJSON
#' @importFrom data.table rbindlist dcast fread setDT fifelse melt data.table as.data.table setnames setorder
#' @importFrom sf st_as_sf
#' @returns `list` containing the following objects:
#' \describe{
#' \item{`data`}{`"tbl_df"` object. Metadata.}
#' \item{`spots`}{`"tbl_df"` object or `"sf"` if `sf == TRUE`. Locations and spot descriptions.}
#' \item{`tags`}{`"tbl_df"` object. Tags and their descriptsions.}
#' \item{`planar`}{Plane elements. Same row IDs as in `data`.}
#' \item{`linear`}{Line elements. Same row IDs as in `data`.}
#' }
#' @name strabo
#' @examples
#' \dontrun{
#' dt <- read_strabo_xls("C:/Users/tstephan/Documents/Lakehead/Field work/StraboSpot_Output_11_01_2022_TS.xlsx")
#' stereoplot()
#' lines(dt$planar, col = "lightgrey", d = 90)
#' points(dt$linear)
#'
#' read_strabo_mobile("C:/Users/tstephan/Documents/Lakehead/Field work/StraboSpot_Search_06_25_2023.txt")
#'
#' read_strabo_JSON("G:/My Drive/Moss_Lake/data.json")
#' }
NULL

#' @rdname strabo
#' @export
read_strabo_xls <- function(file, tag_cols = FALSE, sf = TRUE) {
  # Read data fast with fread (via readxl first, then as data.table)
  data <- readxl::read_xlsx(file, sheet = 1, skip = 2)
  setDT(data)
  setnames(data, make.names(names(data)))

  # Convert Date to POSIXct (quiet, no lubridate)
  data[, Date := as.POSIXct(Date, tz = "UTC")]

  # Calculate Planar.Orientation.Dipdirection (assuming rhr2dd is defined)
  data[, Planar.Orientation.Dipdirection := rhr2dd(Planar.Orientation.Strike)]

  # Linear.Sense logic vectorized and concise
  data[, Linear.Sense := fifelse(Planar.Orientation.Movement == "left_lateral", -1, NA_real_)]
  data[Planar.Orientation.Movement == "right_lateral", Linear.Sense := 1]
  data[Planar.Orientation.Fault.Or.Sz.Type %in% c("sinistral", "reverse"), Linear.Sense := -1]
  data[Planar.Orientation.Fault.Or.Sz.Type %in% c("dextral", "normal", "dextral_normal"), Linear.Sense := 1]

  # If tags columns requested, pivot longer and filter
  if (tag_cols) {
    tag_cols_names <- grep("^Tag", names(data), value = TRUE)
    data <- melt(data, measure.vars = tag_cols_names, variable.name = "Tag", value.name = "temp", variable.factor = FALSE)
    data <- data[temp == "X", !"temp", with = FALSE]
    data[, Tag := sub("^Tag:", "", Tag)]
  }

  # Separate linear and planar data.tables
  data_lines <- data[!is.na(Linear.Orientation.Trend), .SD, .SDcols = patterns("^Linear")]
  data_planes <- data[!is.na(Planar.Orientation.Dipdirection), .SD, .SDcols = patterns("^Planar")]

  # Data without planar or linear cols
  cols_to_exclude <- grep("^(Planar|Linear)", names(data), value = TRUE)
  data0 <- data[, setdiff(names(data), cols_to_exclude), with = FALSE]
  data0[, Planar.Orientation.Unix.Timestamp := data$Planar.Orientation.Unix.Timestamp]

  # Full join on timestamps, using data.table merge with all=TRUE
  res <- merge(data_planes, data_lines,
               by.x = "Planar.Orientation.Unix.Timestamp",
               by.y = "Linear.Orientation.Unix.Timestamp",
               all = TRUE,
               allow.cartesian = TRUE
  )
  res <- merge(res, data0, by = "Planar.Orientation.Unix.Timestamp", all.x = TRUE)

  # Clean up timestamp columns
  res[, Linear.Orientation.Unix.Timestamp := NULL]

  # Construct planes and lines
  planes <- Plane(res$Planar.Orientation.Dipdirection, res$Planar.Orientation.Dip)
  lines <- Line(res$Linear.Orientation.Trend, res$Linear.Orientation.Plunge)

  # Convert to sf if requested
  if (sf) {
    res <- st_as_sf(res, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE, na.fail = FALSE)
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
  # Read the data with data.table::fread for speed
  data0 <- fread(file,
                 sep = "\t", header = TRUE,
                 colClasses = c("integer", rep("character", 3), rep("numeric", 11), "character")
  )

  # Separate lines and planes using data.table syntax
  lines0 <- data0[Type == "L"]
  planes0 <- data0[Type == "P"]

  # Create 'Dipdir' for planes (vectorized)
  planes0[, Dipdir := (`Trd/Strk` + 90) %% 360]

  # Construct line and plane objects from orientation data
  lines <- Line(lines0$`Trd/Strk`, lines0$`Plg/Dip`)
  planes <- Plane(planes0$Dipdir, planes0$`Plg/Dip`)

  # Extract metadata tables, dropping columns as needed
  lines.meta <- lines0[, !c("No.", "Trd/Strk", "Plg/Dip"), with = FALSE]
  planes.meta <- planes0[, !c("No.", "Dipdir", "Trd/Strk", "Plg/Dip"), with = FALSE]

  # Set row names based on 'No.'
  rownames(lines) <- rownames(lines.meta) <- lines0$No.
  rownames(planes) <- rownames(planes.meta) <- planes0$No.

  # Convert metadata to sf objects if requested
  if (sf) {
    lines.meta <- st_as_sf(lines.meta,
                           coords = c("Longitude", "Latitude"),
                           crs = 4326, remove = FALSE, na.fail = FALSE
    )
    planes.meta <- st_as_sf(planes.meta,
                            coords = c("Longitude", "Latitude"),
                            crs = 4326, remove = FALSE, na.fail = FALSE
    )
  }

  # Return list
  list(
    linear = lines,
    linear_data = lines.meta,
    planar = planes,
    planar_data = planes.meta
  )
}


#' @rdname strabo
#' @export
# read_strabo_JSON_old <- function(file, sf = TRUE) {
#   # works but is super slow
#
#   spot.y <- ds_id <- unix_timestamp <- modified_timestamp <- spot.x <- NULL
#
#   time <- tag_id <- tag_name <- NULL
#
#   dat <- rjson::fromJSON(file = file)
#
#   # read tags
#   tags_list <- dat$project$project$tags
#   tags_df <- spot_tags <- tibble()
#   for (t in seq_along(tags_list)) {
#     spots_t <- data.frame(spot = tags_list[[t]]$spots)
#     tags_list[[t]]$spots <- NULL
#     tags_list[[t]]$eon <- NULL
#
#     df_t <- as.data.frame(tags_list[[t]])
#     tags_df <- dplyr::bind_rows(tags_df, df_t)
#     if (nrow(spots_t) > 0) {
#       spots_t$tag <- df_t$id
#       spot_tags <- rbind(spot_tags, spots_t)
#     }
#   }
#   colnames(tags_df) <- paste("tag", colnames(tags_df), sep = "_")
#
#   spot_tags_df <- dplyr::left_join(
#     unique(spot_tags),
#     tags_df |> dplyr::select(tag_id, tag_name),
#     by = c("tag" = "tag_id")
#   ) |>
#     dplyr::mutate(tag_name = paste0("tag:", tag_name)) |>
#     # dplyr::select(-tag) |>
#     tidyr::pivot_wider(names_from = "tag_name", values_from = "tag") |>
#     dplyr::mutate(dplyr::across(!spot, function(x) ifelse(is.na(x), FALSE, TRUE))) |>
#     unique() |>
#     mutate(spot = as.character(spot)) |>
#     rename(spot_id = spot)
#
#   # datasets
#   ds <- seq_along(dat$project$datasets)
#   ds_name <- ds_name_id <- character()
#
#   for (a in ds) {
#     ds_name[a] <- dat$project$datasets[[a]]$name
#     ds_name_id[a] <- names(dat$project$datasets)[[a]]
#   }
#   ds_df <- data.frame(ds_name = ds_name, ds_id = ds_name_id)
#
#   # spot lists
#   spots <- character()
#   for (i in ds_name_id) {
#     ds_spot_ids_i <- dat$project$datasets[[i]]$spotIds
#     spots <- rbind(spots, cbind(ds_id = rep(i, length(ds_spot_ids_i)), spot_id = ds_spot_ids_i))
#   }
#   spots_df <- as.data.frame(spots) |> unique()
#
#   # extract the specified dataset
#   dsn <- data.frame(ds, ds_name) |>
#     # dplyr::filter(ds_name == dataset) |>
#     dplyr::pull(ds)
#   ds_name_id <- names(dat$project$datasets)[dsn]
#
#
#   # ds_spot_ids <- dat$project$datasets[[dsn]]$spotIds # spots in dataset
#   # DSN <- dat$spotsDb[ds_spot_ids]
#   # for (i in seq_along(dat$spotsDb)) {
#   #   dat$spotsDb[[i]]
#   # }
#
#
#   spot_id <- spot <- notes <- date <- character()
#   altitude <- longitude <- latitude <- gps_accuracy <- numeric()
#   orient_df <- data.frame(check.names = FALSE)
#
#   id <- numeric()
#   for (i in seq_along(dat$spotsDb)) {
#     spot_i <- dat$spotsDb[[i]]
#     if (!is.null(spot_i$geometry)) {
#       if (spot_i$geometry$type == "Point" & spot_i$type == "Feature") {
#         spot_id[i] <- spot_i$properties$id
#         spot[i] <- spot_i$properties$name
#         date[i] <- spot_i$properties$date
#         id[i] <- spot_i$properties$id
#
#         notes_i <- spot_i$properties$notes
#         notes[i] <- ifelse(is.null(notes_i), NA, notes_i)
#
#         altitude_i <- spot_i$properties$altitude
#         altitude[i] <- ifelse(is.null(altitude_i), NA, altitude_i)
#
#         longitude[i] <- spot_i$geometry$coordinates[1]
#         latitude[i] <- spot_i$geometry$coordinates[2]
#
#         gps_accuracy_i <- spot_i$properties$gps_accuracy
#         gps_accuracy[i] <- ifelse(is.null(gps_accuracy_i), NA, gps_accuracy_i)
#
#         # orientation data
#         orient_ls_i <- spot_i$properties$orientation_data
#         if (!is.null(orient_ls_i)) {
#           orient_df_i <- tibble()
#           for (j in seq_along(orient_ls_i)) { # loop through each measurement
#             orient_ls_ij <- orient_ls_i[[j]]
#
#             # concatenate directional_indicators into one column
#             if (length(orient_ls_ij$directional_indicators) > 1) {
#               orient_ls_ij$directional_indicators <- toString(orient_ls_ij$directional_indicators)
#             }
#
#             # single P or L measurements
#             if (is.null(orient_ls_ij$associated_orientation)) {
#               # append all measurements of spot
#               # orient_df_i <- plyr::rbind.fill(orient_df_i, as.data.frame(orient_ls_ij, check.names = FALSE)) |>
#               orient_df_i <- plyr::rbind.fill(orient_df_i, as.data.frame(orient_ls_ij, check.names = FALSE))
#               orient_df_i$spot <- spot[i]
#               orient_df_i$spot_id <- spot_id[i]
#               orient_df_i$associated <- FALSE
#             } else {
#               # P + L measurements
#
#               L_ls <- orient_ls_ij$associated_orientation
#               orient_ls_ij$associated_orientation <- NULL
#
#               ## P measurement
#               P <- orient_ls_ij |>
#                 as.data.frame(check.names = FALSE) |>
#                 unique()
#
#               ## L measurement
#               PL <- data.frame(check.names = FALSE)
#               for (k in seq_along(L_ls)) {
#                 if (L_ls[[k]]$type != "planar_orientation") {
#                   L <- L_ls[[k]] |>
#                     as.data.frame(check.names = FALSE) |>
#                     unique()
#                   plunge <- L$plunge
#                   trend <- L$trend
#                   if (!is.na(plunge) & !is.na(trend)) {
#                     L$unix_timestamp <- L$plunge <- L$trend <- NULL # don't need this column twice
#                     colnames(L) <- paste("associated", colnames(L), sep = "_")
#                     L$plunge <- plunge
#                     L$trend <- trend
#
#                     ## combine P + L measurements
#                     PL_k <- cbind(P, L)
#                     PL <- plyr::rbind.fill(PL, PL_k) |>
#                       unique()
#                   }
#                   PL$associated <- TRUE
#                 }
#               }
#
#               # append all P+L measurements of spot
#               if (nrow(orient_df_i) > 0) {
#                 orient_df_i <- plyr::rbind.fill(orient_df_i, PL) |>
#                   unique()
#                 orient_df_i$spot <- spot[i]
#                 orient_df_i$spot_id <- spot_id[i]
#               }
#
#               # orient_ls_ij <- NULL
#             }
#             # combine all spots
#             orient_df <- plyr::rbind.fill(orient_df, orient_df_i)
#           }
#         }
#       }
#     }
#   }
#   fieldbook <- data.frame(spot_id = spot_id, spot = spot, longitude, latitude, altitude, note = notes, time = date, gps_accuracy, id) |>
#     dplyr::mutate(time = lubridate::as_datetime(time)) |>
#     dplyr::arrange(time) |>
#     dplyr::distinct() |>
#     dplyr::left_join(spots_df) |>
#     dplyr::left_join(ds_df) |>
#     dplyr::left_join(spot_tags_df) |>
#     dplyr::select(-c(id))
#   if (sf) fieldbook <- sf::st_as_sf(fieldbook, coords = c("longitude", "latitude"), remove = FALSE, crs = "WGS84", na.fail = FALSE)
#
#   meta <- orient_df |>
#     as.data.frame() |>
#     dplyr::distinct() |>
#     dplyr::left_join(fieldbook, dplyr::join_by("spot_id")) |>
#     dplyr::arrange(time) |>
#     dplyr::select(-c(spot.y, ds_id, spot_id, unix_timestamp, modified_timestamp)) |>
#     dplyr::rename(spot = spot.x)
#   #
#   if (sf) {
#     meta <- sf::st_as_sf(meta)
#   }
#
#   planes <- as.plane(cbind(orient_df$strike + 90, orient_df$dip))
#   lines <- as.line(cbind(orient_df$trend, orient_df$plunge))
#
#   list(
#     data = meta,
#     tags = tags_df,
#     planar = planes,
#     linear = lines
#   )
# }

# read_strabo_JSON_fast <- function(file, sf = TRUE) {
#   # works but it can be faster
#
#   dat <- rjson::fromJSON(file = file)
#
#   ## TAGS
#   tags_list <- dat$project$project$tags %||% list()
#   tags_df <- purrr::map_dfr(tags_list, ~{
#     spots <- tibble::tibble(spot = .x$spots %||% character())
#     id <- .x$id %||% NA
#     name <- .x$name %||% NA
#     if (nrow(spots) > 0) {
#       spots$tag <- id
#       list(tags = dplyr::tibble(tag_id = id, tag_name = name), spot_tags = spots)
#     } else {
#       list(tags = dplyr::tibble(tag_id = id, tag_name = name), spot_tags = NULL)
#     }
#   })
#
#   tags_combined <- bind_rows(purrr::map(tags_list, ~ dplyr::tibble(tag_id = .x$id, tag_name = .x$name)))
#   spot_tags <- bind_rows(purrr::map(tags_list, ~{
#     if (!is.null(.x$spots) && length(.x$spots) > 0) {
#       tibble::tibble(spot = as.character(.x$spots), tag = .x$id)
#     }
#   }))
#
#   if (nrow(spot_tags) > 0) {
#     spot_tags_df <- spot_tags %>%
#       dplyr::left_join(tags_combined, by = c("tag" = "tag_id")) %>%
#       dplyr::mutate(tag_name = paste0("tag:", tag_name)) %>%
#       tidyr::pivot_wider(names_from = "tag_name", values_from = "tag", values_fill = FALSE) %>%
#       dplyr::mutate(across(-spot, ~ !is.na(.x))) %>%
#       dplyr::rename(spot_id = spot)
#   } else {
#     spot_tags_df <- dplyr::tibble(spot_id = character())
#   }
#
#   ## DATASETS
#   datasets <- dat$project$datasets
#   ds_names <- names(datasets)
#   ds_df <- dplyr::tibble(
#     ds_name = purrr::map_chr(datasets, "name", .default = NA_character_),
#     ds_id = ds_names
#   )
#
#   ## SPOTS-DATASET RELATION
#   spots_df <- purrr::imap_dfr(datasets, ~ {
#     dplyr::tibble(
#       ds_id = .y,
#       spot_id = as.character(.x$spotIds %||% character())
#     )
#   }) %>%
#     dplyr::distinct()
#
#   ## SPOTS PROPERTIES
#   spotsDb <- dat$spotsDb %||% list()
#   valid_spots <- purrr::keep(spotsDb, ~ !is.null(.x$geometry) && .x$geometry$type == "Point" && .x$type == "Feature")
#
#   fieldbook_df <- purrr::map_dfr(valid_spots, ~{
#     props <- .x$properties
#     data.frame(
#       spot_id = as.character(props$id %||% NA_character_),   # force character
#       spot = props$name %||% NA_character_,
#       longitude = .x$geometry$coordinates[1] %||% NA_real_,
#       latitude = .x$geometry$coordinates[2] %||% NA_real_,
#       altitude = props$altitude %||% NA_real_,
#       note = props$notes %||% NA_character_,
#       time = props$date %||% NA_character_,
#       gps_accuracy = props$gps_accuracy %||% NA_real_,
#       id = as.character(props$id %||% NA_character_)         # consistent type
#     )
#   }) %>%
#     dplyr::mutate(time = as.POSIXct(time, tz = "UTC", tryFormats = c("%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%d %H:%M:%S"))) |>
#     dplyr::left_join(spots_df, by = "spot_id") %>%
#     dplyr::left_join(ds_df, by = "ds_id") %>%
#     dplyr::left_join(spot_tags_df, by = "spot_id") %>%
#     dplyr::arrange(time) %>%
#     dplyr::distinct() %>%
#     dplyr::select(-id)
#
#   if (sf) {
#     fieldbook_df <- sf::st_as_sf(fieldbook_df, coords = c("longitude", "latitude"), crs = "WGS84", remove = FALSE)
#   }
#
#   ## ORIENTATION DATA (extract only if needed)
#   orient_df <- purrr::map_dfr(valid_spots, ~ {
#     orient <- .x$properties$orientation_data
#     if (!is.null(orient)) {
#       purrr::map_dfr(orient, function(o) {
#         if (!is.null(o$associated_orientation)) {
#           L <- purrr::map_dfr(o$associated_orientation, function(l) {
#             if (l$type != "planar_orientation") {
#               dplyr::tibble(
#                 associated_trend = l$trend %||% NA_real_,
#                 associated_plunge = l$plunge %||% NA_real_
#               )
#             }
#           })
#           P <- dplyr::tibble(
#             strike = o$strike %||% NA_real_,
#             dip = o$dip %||% NA_real_
#           )
#           bind_cols(P, L)
#         } else {
#           dplyr::tibble(
#             strike = o$strike %||% NA_real_,
#             dip = o$dip %||% NA_real_
#           )
#         }
#       })
#     }
#   })
#
#   if (nrow(orient_df) > 0) {
#     # if (sf) orient_df <- sf::st_as_sf(orient_df)
#     planes <- as.plane(cbind(orient_df$strike + 90, orient_df$dip))
#     lines <- as.line(cbind(orient_df$associated_trend, orient_df$associated_plunge))
#   } else {
#     planes <- NULL
#     lines <- NULL
#   }
#
#   list(
#     data = fieldbook_df,
#     tags = tags_combined,
#     planar = planes,
#     linear = lines
#   )
# }

read_strabo_JSON <- function(file, sf = TRUE) {
  # --- Load JSON ---
  dat <- rjson::fromJSON(file = file)

  # --- Tags ---
  tags_list <- dat$project$project$tags
  tags_dt <- rbindlist(
    lapply(tags_list, function(tag) {
      tag$eon <- NULL
      spots <- data.table(spot_id = unlist(tag$spots))
      tag$spots <- NULL
      tag_dt <- as.data.table(tag)
      if (nrow(spots) > 0) {
        spots[, tag_id := tag_dt$id]
        return(spots)
      }
      NULL
    }),
    fill = TRUE
  )

  tag_info_dt <- rbindlist(
    lapply(tags_list, function(tag) {
      tag$spots <- NULL
      tag$eon <- NULL
      as.data.table(tag)
    }),
    fill = TRUE
  )
  setnames(tag_info_dt, paste0("tag_", names(tag_info_dt)))

  spot_tags_dt <- merge(
    unique(tags_dt),
    tag_info_dt[, .(tag_id = tag_id, tag_name = tag_name)],
    by.x = "tag_id", by.y = "tag_id",
    all.x = TRUE
  )
  spot_tags_dt[, tag_col := paste0("tag:", tag_name)]
  spot_tags_wide <- dcast(spot_tags_dt[, .(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)],
                          spot_id ~ tag_name,
                          fill = FALSE
  )

  # --- Datasets ---
  datasets <- dat$project$datasets
  ds_dt <- data.table(
    ds_name = sapply(datasets, `[[`, "name"),
    ds_id = names(datasets)
  )

  # Spot-dataset mapping
  spots_dt <- rbindlist(
    lapply(names(datasets), function(ds_id) {
      spot_ids <- datasets[[ds_id]]$spotIds
      if (is.null(spot_ids) || length(spot_ids) == 0) {
        # Return empty data.table with correct column types
        return(data.table(ds_id = character(0), spot_id = character(0)))
      }
      data.table(ds_id = ds_id, spot_id = spot_ids)
    }),
    fill = TRUE
  )
  spots_dt <- unique(spots_dt)

  # --- Spots ---
  spots_db <- dat$spotsDb

  fieldbook_list <- lapply(spots_db, function(spot) {
    if (!is.null(spot$geometry) && spot$geometry$type == "Point" && spot$type == "Feature") {
      props <- spot$properties
      data.table(
        spot_id = props$id %||% NA_character_,
        spot = props$name %||% NA_character_,
        longitude = spot$geometry$coordinates[1] %||% NA_real_,
        latitude = spot$geometry$coordinates[2] %||% NA_real_,
        altitude = props$altitude %||% NA_real_,
        note = props$notes %||% NA_character_,
        time = as.POSIXct(props$date %||% NA_character_, tz = "UTC", tryFormats = c("%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%d %H:%M:%S")),
        gps_accuracy = props$gps_accuracy %||% NA_real_,
        id = props$id %||% NA_character_
      )
    } else {
      NULL
    }
  })
  fieldbook_dt <- rbindlist(fieldbook_list, fill = TRUE)

  fieldbook_dt[, spot_id := as.character(spot_id)]
  spots_dt[, spot_id := as.character(spot_id)]
  spot_tags_wide[, spot_id := as.character(spot_id)]


  # --- Join with dataset and tag info ---
  fieldbook_dt <- merge(fieldbook_dt, spots_dt, by = "spot_id", all.x = TRUE)
  fieldbook_dt <- merge(fieldbook_dt, ds_dt, by = "ds_id", all.x = TRUE)
  fieldbook_dt <- merge(fieldbook_dt, spot_tags_wide, by = "spot_id", all.x = TRUE)
  fieldbook_dt <- unique(fieldbook_dt)
  setorder(fieldbook_dt, time)
  fieldbook_dt[, id := NULL] # remove id column

  if (sf) {
    fieldbook_dt <- st_as_sf(fieldbook_dt, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  }

  # --- Orientation data ---
  orient_list <- lapply(spots_db, function(spot) {
    if (!is.null(spot$geometry) && spot$geometry$type == "Point" && spot$type == "Feature") {
      props <- spot$properties
      spot_id <- props$id
      spot_name <- props$name
      orient_data <- props$orientation_data
      if (!is.null(orient_data)) {
        rbindlist(lapply(orient_data, function(orient) {
          assoc <- orient$associated_orientation
          orient$associated_orientation <- NULL
          orient_dt <- as.data.table(orient)
          orient_dt[, `:=`(spot = spot_name, spot_id = spot_id, associated = FALSE)]
          if (!is.null(assoc)) {
            assoc_dt <- rbindlist(lapply(assoc, function(a) {
              if (a$type != "planar_orientation") {
                a_dt <- as.data.table(a)
                setnames(a_dt, names(a_dt), paste0("associated_", names(a_dt)))
                return(a_dt)
              }
              NULL
            }), fill = TRUE)
            if (nrow(assoc_dt) > 0) {
              cbind(orient_dt, assoc_dt, associated = TRUE)
            } else {
              orient_dt
            }
          } else {
            orient_dt
          }
        }), fill = TRUE)
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  orient_dt <- rbindlist(orient_list, fill = TRUE)
  orient_dt <- unique(orient_dt)

  if (nrow(orient_dt) > 0) {
    # if (sf) orient_df <- sf::st_as_sf(orient_df)
    planes <- Plane(orient_dt$strike + 90, orient_dt$dip)
    lines <- Line(orient_dt$associated_trend, orient_dt$associated_plunge)
  } else {
    planes <- NULL
    lines <- NULL
  }


  # --- Return ---
  list(
    data = orient_dt,
    spots = fieldbook_dt,
    tags = tag_info_dt,
    planar = planes,
    linear = lines
  )
}
