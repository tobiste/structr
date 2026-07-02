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
#' 
#' @importFrom readxl read_xlsx
#' @importFrom rjson fromJSON
# #' @importFrom data.table rbindlist dcast fread setDT fifelse melt data.table as.data.table setnames setorder setcolorder patterns .SD .( :=
#' @import data.table
#' @importFrom sf st_as_sf
#' 
#' @returns `list` containing the following objects:
#' \describe{
#' \item{`data`}{`"tbl_df"` object. Metadata.}
#' \item{`spots`}{`"tbl_df"` object or `"sf"` if `sf == TRUE`. Locations and spot descriptions.}
#' \item{`tags`}{`"tbl_df"` object. Tags and their descriptions.}
#' \item{`planar`}{Plane elements. Same row IDs as in `data`.}
#' \item{`linear`}{Line elements. Same row IDs as in `data`.}
#' }
#' @name strabo
#'
#' @family strabo
#' @seealso [drillcore_transformation()]
#'
#' @examples
#' \dontrun{
#' # import from excel file
#' read_strabo_xlsx("path/to/my/file.xlsx")
#'
#' # import from text file
#' read_strabo_mobile("path/to/my/file.txt")
#'
#' # import from .json file
#' read_strabo_JSON("path/to/my/file.json")
#' }
NULL

#' @rdname strabo
#' @export
read_strabo_xlsx <- function(file, tag_cols = FALSE, sf = TRUE) {
  Date <- Dipdir <- Linear.Orientation.Trend <- Linear.Orientation.Unix.Timestamp <- Linear.Sense <- Planar.Orientation.Dipdirection <-
    temp <- Planar.Orientation.Fault.Or.Sz.Type <- Planar.Orientation.Movement <- Planar.Orientation.Strike <- Planar.Orientation.Unix.Timestamp <-
    Linear.Orientation.Plunge <- Tabular.Orientation.Dip <- Orientation.Unix.Timestamp <-
    Tag <- NULL


  # Read data fast with fread (via readxl first, then as data.table)
  data <- readxl::read_xlsx(file, sheet = 1, skip = 2)
  setDT(data)
  setnames(data, make.names(names(data)))

  # Convert Date to POSIXct (quiet, no lubridate)
  data[, Date := as.POSIXct(Date, tz = "UTC")]


  # old_names <- names(data)[startsWith(names(data), "Tabular.")]
  # setnames(data, old = old_names, new = gsub("^Tabular.", "Planar.", old_names))


  # Linear.Sense logic vectorized and concise
  data[, Linear.Sense := fifelse(Planar.Orientation.Movement == "left_lateral", -1, NA_real_)]
  data[Planar.Orientation.Movement == "right_lateral", Linear.Sense := 1]
  data[Planar.Orientation.Fault.Or.Sz.Type %in% c("sinistral", "reverse"), Linear.Sense := -1]
  data[Planar.Orientation.Fault.Or.Sz.Type %in% c("dextral", "normal", "dextral_normal"), Linear.Sense := 1]


  if (!"Planar.Orientation.Unix.Timestamp" %in% names(data) &&
    "Planar.Orientation.Modified.Timestamp" %in% names(data)) {
    setnames(data, "Planar.Orientation.Modified.Timestamp", "Planar.Orientation.Unix.Timestamp")
  }
  if (!"Linear.Orientation.Unix.Timestamp" %in% names(data) &&
    "Linear.Orientation.Modified.Timestamp" %in% names(data)) {
    setnames(data, "Linear.Orientation.Modified.Timestamp", "Linear.Orientation.Unix.Timestamp")
  }
  if (!"Tabular.Orientation.Unix.Timestamp" %in% names(data) &&
    "Tabular.Orientation.Modified.Timestamp" %in% names(data)) {
    setnames(data, "Tabular.Orientation.Modified.Timestamp", "Tabular.Orientation.Unix.Timestamp")
  }

  # If tags columns requested, pivot longer and filter
  if (tag_cols) {
    tag_cols_names <- grep("^Tag", names(data), value = TRUE)
    data <- melt(data, measure.vars = tag_cols_names, variable.name = "Tag", value.name = "temp", variable.factor = FALSE)
    data <- data[temp == "X", !"temp", with = FALSE]
    data[, Tag := sub("^Tag:", "", Tag)]
  }

  # Separate linear and planar data.tables
  data_lines <- data[!is.na(Linear.Orientation.Plunge), .SD, .SDcols = data.table::patterns("^Linear")]
  data_planes <- data[!is.na(Planar.Orientation.Dipdirection), .SD, .SDcols = data.table::patterns("^Planar")]
  data_tabular <- data[!is.na(Tabular.Orientation.Dip), .SD, .SDcols = data.table::patterns("^Tabular")]
  old_names <- grep("^Tabular\\.", names(data_tabular), value = TRUE)
  setnames(data_tabular, old = old_names, new = gsub("^Tabular\\.", "Planar.", old_names))
  data_planes <- rbindlist(list(data_planes, data_tabular), fill = TRUE)

  # Data without planar or linear cols
  cols_to_exclude <- grep("^(Planar|Linear|Tabular)", names(data), value = TRUE)
  data0 <- data[, setdiff(names(data), cols_to_exclude), with = FALSE]

  t <- ifelse(is.na(data$Linear.Orientation.Unix.Timestamp), data$Planar.Orientation.Unix.Timestamp, data$Linear.Orientation.Unix.Timestamp)
  t <- ifelse(is.na(t), data$Tabular.Orientation.Unix.Timestamp, t)

  data0[, Orientation.Unix.Timestamp := t]
  data0 <- data0[!is.na(Orientation.Unix.Timestamp), ]

  # Full join on timestamps, using data.table merge with all=TRUE
  res <- merge(data_planes, data_lines,
    by.x = "Planar.Orientation.Unix.Timestamp",
    by.y = "Linear.Orientation.Unix.Timestamp",
    all = TRUE,
    allow.cartesian = TRUE
  )
  res <- merge(res, data0, by.x = "Planar.Orientation.Unix.Timestamp", by.y = "Orientation.Unix.Timestamp", all.x = TRUE)
  setnames(res, "Planar.Orientation.Unix.Timestamp", "Orientation.Unix.Timestamp")

  # res2 <- merge(data_tabular, data_lines,
  #              by.x = "Tabular.Orientation.Unix.Timestamp",
  #              by.y = "Linear.Orientation.Unix.Timestamp",
  #              all = TRUE,
  #              allow.cartesian = TRUE
  # )
  # res2 <- merge(res2, data0, by.x = "Tabular.Orientation.Unix.Timestamp", by.y = "Orientation.Unix.Timestamp", all.x = TRUE)
  # setnames(res2, "Tabular.Orientation.Unix.Timestamp", "Orientation.Unix.Timestamp")
  #
  # res <- rbindlist(list(res1, res2), fill = TRUE)

  # Clean up timestamp columns
  # if(!is.null(res$Linear.Orientation.Unix.Timestamp)) res[, Linear.Orientation.Unix.Timestamp := NULL]

  # Calculate Planar.Orientation.Dipdirection (assuming rhr2dd is defined)
  # data[, Planar.Orientation.Dipdirection := rhr2dd(Planar.Orientation.Strike)]

  # Construct planes and lines
  planes <- Plane(rhr2dd(res$Planar.Orientation.Strike), res$Planar.Orientation.Dip)
  lines <- Line(res$Linear.Orientation.Trend, res$Linear.Orientation.Plunge)

  # Convert to sf if requested
  if (sf) {
    res <- st_as_sf(res, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE, na.fail = FALSE)
  }

  ls <- list(
    data = res,
    planar = planes,
    linear = lines
  )
  as.strabo(ls)
}

#' @rdname strabo
#' @export
read_strabo_xls <- read_strabo_xlsx

#' @rdname strabo
#' @export
read_strabo_mobile <- function(file, sf = TRUE) {
  `Plg/Dip` <- `Trd/Strk` <- Dipdir <- Type <- Strk <- NULL
  # Read the data with data.table::fread for speed
  data0 <- fread(file,
    sep = "\t", header = TRUE,
    colClasses = c("integer", rep("character", 3), rep("numeric", 11), "character")
  )

  data <- data0
  data$`Trd/Strk` <- NULL
  data$`Plg/Dip` <- NULL


  data$azimuth <- ifelse(data0$Type == "L", data0$`Trd/Strk`, NA_real_)
  strike <- ifelse(data0$Type == "P", data0$`Trd/Strk`, NA_real_)
  data$plunge <- ifelse(data0$Type == "L", data0$`Plg/Dip`, NA_real_)
  data$dip <- ifelse(data0$Type == "P", data0$`Plg/Dip`, NA_real_)
  data$dip_direction <- rhr2dd(strike)

  # Separate lines and planes using data.table syntax
  # lines0 <- data0[Type == "L"]
  # planes0 <- data0[Type == "P"]

  # Create 'Dipdir' for planes (vectorized)
  # planes0[, Dipdir := (`Trd/Strk` + 90) %% 360]

  # Construct line and plane objects from orientation data
  # lines <- Line(lines0$`Trd/Strk`, lines0$`Plg/Dip`)
  # planes <- Plane(planes0$Dipdir, planes0$`Plg/Dip`)
  lines <- Line(data$azimuth, data$plunge)
  planes <- Plane(data$dip_direction, data$dip)

  # Extract metadata tables, dropping columns as needed
  # lines.meta <- lines0[, !c("No.", "Trd/Strk", "Plg/Dip"), with = FALSE]
  # planes.meta <- planes0[, !c("No.", "Dipdir", "Trd/Strk", "Plg/Dip"), with = FALSE]

  # Set row names based on 'No.'
  # rownames(lines) <- rownames(lines.meta) <- lines0$No.
  # rownames(planes) <- rownames(planes.meta) <- planes0$No.

  # Convert metadata to sf objects if requested
  if (sf) {
    # lines.meta <- st_as_sf(lines.meta,
    #   coords = c("Longitude", "Latitude"),
    #   crs = 4326, remove = FALSE, na.fail = FALSE
    # )
    # planes.meta <- st_as_sf(planes.meta,
    #   coords = c("Longitude", "Latitude"),
    #   crs = 4326, remove = FALSE, na.fail = FALSE
    # )
    data <- st_as_sf(data,
      coords = c("Longitude", "Latitude"),
      crs = 4326, remove = FALSE, na.fail = FALSE
    )
  }

  # Return list
  ls <- list(
    data = data,
    # linear_data = lines.meta,
    planar = planes,
    linear = lines
    # planar_data = planes.meta
  )
  as.strabo(ls)
}


#' @rdname strabo
#' @export
read_strabo_JSON <- function(file, sf = TRUE) {
  tag_name <- spot_id <- tag_col <- tag_id <- time <- id <- NULL
  associated <- type <- associated_type <- unix_timestamp <- NULL

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
  tags_dt$spot_id <- as.character(tags_dt$spot_id)

  if (is.null(tags_list)) {
    tag_info_dt <- NULL
  } else {
    tag_info_dt <- rbindlist(
      lapply(tags_list, function(tag) {
        tag$spots <- NULL
        tag$eon <- NULL
        as.data.table(tag)
      }),
      fill = TRUE
    )
    setnames(tag_info_dt, paste0("tag_", names(tag_info_dt)))
  }

  # spot_tags_dt <- merge(
  #   unique(tags_dt),
  #   tag_info_dt[, .(tag_id = tag_id, tag_name = tag_name)],
  #   by.x = "tag_id", by.y = "tag_id",
  #   all.x = TRUE
  # )

  if (nrow(unique(tags_dt)) == 0) {
    spot_tags_dt <- tag_info_dt[, list(tag_id = tag_id, tag_name = tag_name)]
  } else {
    spot_tags_dt <- merge(
      unique(tags_dt),
      tag_info_dt[, list(tag_id = tag_id, tag_name = tag_name)],
      by.x = "tag_id", by.y = "tag_id",
      all.x = TRUE
    )
  }

  if (!is.null(tags_list)) spot_tags_dt[, tag_col := paste0("tag:", tag_name)]

  # spot_tags_wide <- dcast(spot_tags_dt[, .(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)],
  #   spot_id ~ tag_name,
  #   fill = FALSE
  # )
  if (is.null(spot_tags_dt$spot_id)) {
    spot_tags_wide <- spot_tags_dt[, list(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)]
    spot_tags_wide$spot_id <- NA_character_
  } else {
    spot_tags_wide <- dcast(
      spot_tags_dt[, list(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)],
      spot_id ~ tag_name,
      fill = FALSE
    )
  }

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
  spots_dt$spot_id <- as.character(spots_dt$spot_id)

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
  if (!is.null(tags_list)) spot_tags_wide[, spot_id := as.character(spot_id)]


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
                # a_dt$associated <- NULL # drop associated column in lineations
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
  orient_dt$spot_id <- as.character(orient_dt$spot_id)

  # Fix duplicate associated columns
  idx <- which(names(orient_dt) == "associated")
  # orient_dt[, associated :=
  #      Reduce(`|`, lapply(.SD, coalesce, FALSE)),
  #    .SDcols = idx]
  orient_dt[, associated := Reduce(`|`, lapply(idx, function(i) {
    col <- orient_dt[[i]]
    replace(col, is.na(col), FALSE)
  }))]
  # Remove the duplicate column (now at its original position)
  orient_dt[, (idx[2]) := NULL]

  # add strabo measurement column numbers if not present:
  expected_cols <- c(
    "associated", "associated_defined_by", "associated_feature_type", "associated_id",
    "associated_label", "associated_notes", "associated_plunge", "associated_quality",
    "associated_trend", "associated_type", "defined_by", "dip",
    "dip_direction", "directional_indicators", "fault_or_sz_type", "feature_type",
    "foliation_defined_by", "foliation_type", "id", "label",
    "length", "modified_timestamp", "movement", "movement_amount_m",
    "movement_amount_qualifier", "movement_justification", "notes", "other_movement",
    "other_movement_justification", "plunge", "quality", "rake",
    "spot", "spot_id", "strike", "tabularity",
    "thickness", "trend", "type", "unix_timestamp",
    "vein_fill", "vein_type"
  )
  missing_cols <- setdiff(expected_cols, names(orient_dt))
  orient_dt[, (missing_cols) := NA]

  orient_dt[, associated := fifelse(is.na(orient_dt$associated), FALSE, associated)]
  # orient_dt$associated <- ifelse(is.na(orient_dt$associated), FALSE, orient_dt$associated)

  # select only planes that are not associated with a line
  only_planes <- orient_dt[type == "planar_orientation" | type == "tabular_orientation" & !associated, ]
  only_planes[, c("trend", "plunge", grep("^associated_", names(only_planes), value = TRUE)) := NULL] # drop all line and associated related columns
  setnames(only_planes, old = "type", new = "planar_type") # rename type columns
  only_planes[, which(sapply(only_planes, \(x) all(is.na(x)))) := NULL] # drop all empty columns

  # select all lines that are not associated with a plane
  only_lines <- orient_dt[type == "linear_orientation" & !associated, ]
  only_lines[, grep("^associated_", names(only_lines), value = TRUE) := NULL] # drop all associated related columns
  line_specifiers1 <- c("defined_by", "feature_type", "label", "notes", "quality", "type")
  setnames(only_lines, old = line_specifiers1, new = paste0("linear_", line_specifiers1)) # rename all line related columns
  only_lines[, which(sapply(only_lines, \(x) all(is.na(x)))) := NULL] # drop all empty columns

  # select all associated planes and lines
  only_associated_lines <- orient_dt[associated_type == "linear_orientation" & associated, ]
  drop_cols1 <- c("trend", "plunge")
  only_associated_lines[, (drop_cols1) := NULL] # drop trend and plunge columns
  old_names <- names(only_associated_lines)[startsWith(names(only_associated_lines), "associated_")]
  setnames(only_associated_lines, old = old_names, new = gsub("^associated_", "linear_", old_names))
  setnames(only_associated_lines, old = c("linear_trend", "linear_plunge"), new = c("trend", "plunge"))
  only_associated_lines[, which(sapply(only_associated_lines, \(x) all(is.na(x)))) := NULL] # drop all empty columns


  orient_dt2 <- rbindlist(list(only_planes, only_lines, only_associated_lines), fill = TRUE)
  
  
  if("strike" %in% names(orient_dt2)){
    orient_dt2$dip_direction <- rhr2dd(orient_dt2$strike)
  } else if("dip_direction" %in% names(orient_dt2)){
    orient_dt2$strike <- dd2rhr(orient_dt2$dip_direction)
  }
  
  setcolorder(orient_dt2, c("id", "dip_direction", "dip", "strike", "trend", "plunge", "associated"))

  # drop all empty columns
  orient_dt2[, which(sapply(orient_dt2, \(x) all(is.na(x)))) := NULL]

  # sort the data
  setorder(orient_dt2, unix_timestamp)

  # remove duplicated planes
  dup_ids <- orient_dt2[duplicated(id) | duplicated(id, fromLast = TRUE), unique(id)]
  orient_dt2 <- orient_dt2[!(id %in% dup_ids) | (!is.na(dip) & !is.na(plunge))]
  orient_dt2 <- unique(orient_dt2, by = c("id", "dip", "dip_direction"))
  orient_dt2 <- unique(orient_dt2, by = c("id", "plunge", "trend"))
  
  
  if (nrow(orient_dt2) > 0) {
    planes <- Plane(rhr2dd(orient_dt2$strike), orient_dt2$dip)
    lines <- Line(orient_dt2$trend, orient_dt2$plunge)
  } else {
    planes <- NULL
    lines <- NULL
  }


  # --- Return ---
  ls <- list(
    data = orient_dt2,
    spots = fieldbook_dt,
    tags = tag_info_dt,
    planar = planes,
    linear = lines
  )
  strabo <- as.strabo(ls)
  structure(strabo, class = append(class(strabo), "strabo.json"))
}

# if (getRversion() >= "2.15.1")  utils::globalVariables(".")

is.strabo <- function(x) inherits(x, "strabo")

as.strabo <- function(x) {
  structure(x, class = append(class(x), "strabo"))
}

# #' Print strabo imports
# #'
# #' @param n integer. Number of rows to show for each data set.
# #'
# #' @returns combined objects as class defined in `.class`
# #' @method print strabo
# #' @exportS3Method base::print
# #'
# #' @importFrom utils head
# #'
# #' @keywords internal
# #'
# #' @noRd
# print.strabo <- function(x, ..., n = 6L){
#   list(
#     data = utils::head(x$data, n = n, ...),
#     spots = utils::head(x$spots, n = n, ...),
#     tags = utils::head(x$tags, n = n, ...),
#     planar = x$planar,
#     linear = x$linear
#   )
# }


#' Subsetting StraboSpot Projects
#' 
#' Returns a subset from `"strabo"` objects which meet conditions.
#' 
#' @param x object of class `"strabo"`
#' @param ds character. Dataset that should be the base for subsetting. 
#'  One of `"data"` (the default), `"spots"`, or `"tags"`.
#' @param ... arguments to be passed to [subset()]
#' 
#' @returns An object of class `"strabo"` containing just the selected rows and 
#' columns of `"x"`.
#' 
#' @exportS3Method base::subset
#' 
#' @family strabo
#' 
#' @seealso [subset()]
#' 
#' @examples
#' subset(strabo_prj, strabo_prj$data$quality > 3)
subset.strabo <- function(x, ..., ds = c("data", "spots", "tags")) {
  newrowid <- spot_id <- NULL
  
  if (inherits(x, "strabo.json")) {
    ds <- match.arg(ds)
    if (ds == "spots") {
      x$spots$newrowid <- seq_len(nrow(x$spots))
      x_subset <- subset(x$spots, ...)
      selected_rows <- x_subset$newrowid
      x_subset$newrowid <- NULL

      ls <- list(
        data = subset(x$data, ...),
        spots = x_subset,
        planar = x$planar[selected_rows, ],
        linear = x$linear[selected_rows, ]
      )
    } else if (ds == "tags") {
      tags_data_merged <- merge(x$tags, x$data, by.x = "spots_id", by.y = "spot_id", all.x = TRUE)
      tags_subset <- subset(x$tags, ...)
      merged_subset <- subset(tags_data_merged, ...)


      x$data$newrowid <- seq_len(nrow(x$data))
      x_subset <- subset(x$data, spot_id %in% merged_subset$spot_id)
      selected_rows <- x_subset$newrowid
      x_subset$newrowid <- NULL

      ls <- list(
        data = x_subset,
        spots = subset(x$spots, spot_id %in% merged_subset$spot_id),
        tags = tags_subset,
        planar = x$planar[selected_rows, ],
        linear = x$linear[selected_rows, ]
      )
    } else { # ds == 'data'
      x$data$newrowid <- seq_len(nrow(x$data))
      x_subset <- subset(x$data, ...)
      selected_rows <- x_subset$newrowid
      x_subset$newrowid <- NULL

      ls <- list(
        data = x_subset,
        spots = subset(x$spots, spot_id %in% x_subset$spot_id),
        tags = x$tags,
        planar = x$planar[selected_rows, ],
        linear = x$linear[selected_rows, ]
      )
    }
  } else {
    x$data$newrowid <- seq_len(nrow(x$data))
    x_subset <- subset(x$data, ...)
    selected_rows <- x_subset$newrowid
    x_subset$newrowid <- NULL

    ls <- list(
      data = x_subset,
      planar = x$planar[selected_rows, ],
      linear = x$linear[selected_rows, ]
    )
  }

  as.strabo(ls)
}

#' Sorting StraboSpot Projects
#' 
#' Sorts `strabo` objects in the order determined by one or more other objects, typically columns in the
#' `data` data.frame. This will also sort the `planar` and `linear` objects.
#' 
#' @param x object of class `"strabo"`
#' @param y Variables to sort by. Typically a column in `x$data` data.frame. To
#' sort by more than one column, take the form `list(g1, g2, ...)`.
#' @param ... Additional arguments, typically passed on to [order()].
#'  
#' @returns A sorted version of `x`. 
#'
#' @exportS3Method base::sort_by
#' 
#' @examples
#' sort_by(strabo_prj, strabo_prj$data$dip_direction)
sort_by.strabo <- function(x, y, ...){
  # spot_id <- NULL
  newrowid <- NULL
  
  
  x$data$newrowid <- seq_len(nrow(x$data))
  x_sort <- sort_by(x$data, y, ...)
  selected_rows <- x_sort$newrowid
  x_sort$newrowid <- NULL
  
  if (inherits(x, "strabo.json")) {
  ls <- list(
    data = x_sort,
    spots = x$spots,
    tags = x$tags,
    planar = x$planar[selected_rows, ],
    linear = x$linear[selected_rows, ]
  )
  } else {
    ls <- list(
      data = x_sort,
      planar = x$planar[selected_rows, ],
      linear = x$linear[selected_rows, ]
    )
  }
  
  as.strabo(ls)
}
