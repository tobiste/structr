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
# #' @importFrom data.table rbindlist dcast fread setDT fifelse melt data.table as.data.table setnames setorder patterns .SD .( :=
#' @import data.table
#' @importFrom sf st_as_sf
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
#' @seealso [drillcore_transformation()]
#'
#' @examples
#' \dontrun{
#' # import from excel file
#' read_strabo_xls("path/to/my/file.xlsx")
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
read_strabo_xls <- function(file, tag_cols = FALSE, sf = TRUE) {
  Date <- Dipdir <- Linear.Orientation.Trend <- Linear.Orientation.Unix.Timestamp <- Linear.Sense <- Planar.Orientation.Dipdirection <-
    temp <- Planar.Orientation.Fault.Or.Sz.Type <- Planar.Orientation.Movement <- Planar.Orientation.Strike <- Planar.Orientation.Unix.Timestamp <- Tag <- NULL


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
  data_lines <- data[!is.na(Linear.Orientation.Trend), .SD, .SDcols = data.table::patterns("^Linear")]
  data_planes <- data[!is.na(Planar.Orientation.Dipdirection), .SD, .SDcols = data.table::patterns("^Planar")]

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
  if(!is.null(res$Linear.Orientation.Unix.Timestamp)) res[, Linear.Orientation.Unix.Timestamp := NULL]

  # Construct planes and lines
  planes <- Plane(res$Planar.Orientation.Dipdirection, res$Planar.Orientation.Dip)
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
  #lines0 <- data0[Type == "L"]
  #planes0 <- data0[Type == "P"]

  # Create 'Dipdir' for planes (vectorized)
  #planes0[, Dipdir := (`Trd/Strk` + 90) %% 360]

  # Construct line and plane objects from orientation data
  #lines <- Line(lines0$`Trd/Strk`, lines0$`Plg/Dip`)
  #planes <- Plane(planes0$Dipdir, planes0$`Plg/Dip`)
  lines <- Line(data$azimuth, data$plunge)
  planes <- Plane(data$dip_direction, data$dip)
  
  # Extract metadata tables, dropping columns as needed
  #lines.meta <- lines0[, !c("No.", "Trd/Strk", "Plg/Dip"), with = FALSE]
  #planes.meta <- planes0[, !c("No.", "Dipdir", "Trd/Strk", "Plg/Dip"), with = FALSE]

  # Set row names based on 'No.'
  #rownames(lines) <- rownames(lines.meta) <- lines0$No.
  #rownames(planes) <- rownames(planes.meta) <- planes0$No.

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

  if(is.null(tags_list)){
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
  
  if(nrow(unique(tags_dt)) == 0){
    spot_tags_dt <- tag_info_dt[, list(tag_id = tag_id, tag_name = tag_name)]
  } else {
  spot_tags_dt <- merge(
    unique(tags_dt),
    tag_info_dt[, list(tag_id = tag_id, tag_name = tag_name)],
    by.x = "tag_id", by.y = "tag_id",
    all.x = TRUE
  )
  }

  if(!is.null(tags_list)) spot_tags_dt[, tag_col := paste0("tag:", tag_name)]

  # spot_tags_wide <- dcast(spot_tags_dt[, .(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)],
  #   spot_id ~ tag_name,
  #   fill = FALSE
  # )
  if(is.null(spot_tags_dt$spot_id)){
    spot_tags_wide <- spot_tags_dt[, list(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)]
    spot_tags_wide$spot_id <- NA_character_
  } else {
  spot_tags_wide <- dcast(
    spot_tags_dt[, list(spot_id, tag_name = paste0("tag:", tag_name), value = TRUE)],
    spot_id ~ tag_name,
    fill = FALSE
  )}

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
  if(!is.null(tags_list)) spot_tags_wide[, spot_id := as.character(spot_id)]


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
  orient_dt[, associated :=
       Reduce(`|`, lapply(.SD, coalesce, FALSE)),
     .SDcols = idx]
  orient_dt[, idx[-1] := NULL] # remove second associate column
  

  
  if(any(orient_dt$associated)){
  # Condition: trend is NA but associated_trend is not
  fill_from_associated <- is.na(orient_dt$trend) & !is.na(orient_dt$associated_trend)
  
  orient_dt[fill_from_associated, `:=`(
    plunge       = associated_plunge,
    trend        = associated_trend,
    linear_type  = associated_feature_type
  )]
  
  # linear_type for linear_orientation rows (where not already set from above)
  orient_dt[!fill_from_associated & type == "linear_orientation",
            linear_type := feature_type]
  
  # planar_type for planar_orientation rows
  orient_dt[type == "planar_orientation",
            planar_type := feature_type]
  
  # Drop columns
  drop_cols <- c(
    "type", "feature_type", "associated_unix_timestamp",
    "associated_feature_type", "associated_id", "associated_type",
    "associated_trend", "associated_plunge"
  )
  } else {
    orient_dt[,linear_type := feature_type]
    
    # planar_type for planar_orientation rows
    orient_dt[type == "planar_orientation",
              planar_type := feature_type]
    
    # Drop columns
    drop_cols <- c(
      "type", "feature_type"
    )
  }
  orient_dt[, (drop_cols) := NULL]
  
  
  if (nrow(orient_dt) > 0) {
    planes <- Plane(rhr2dd(orient_dt$strike), orient_dt$dip)
    lines <- Line(orient_dt$trend, orient_dt$plunge)
  } else {
    planes <- NULL
    lines <- NULL
  }


  # --- Return ---
  ls <- list(
    data = orient_dt,
    spots = fieldbook_dt,
    tags = tag_info_dt,
    planar = planes,
    linear = lines
  )
  as.strabo(ls)
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

#' @exportS3Method base::subset
subset.strabo <- function(x, ...){
  #stopifnot(is.structr(x))
  x$data$newrowid <- seq_len(nrow(x$data))
  x_subset <- subset(x$data, ...)
  selected_rows <- x_subset$newrowid
  x_subset$newrowid <- NULL
  
  list(
    data = x_subset,
    planar = x$planar[selected_rows, ],
    linear = x$linear[selected_rows, ]
  ) |> 
    as.strabo()
}