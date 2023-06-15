#' Import data from StraboSpot export
#' 
#' @param file the name of the file which the data are to be read from. 
#' @param tag_cols logical. Whether the Tag columns should be summarized in a 
#' single column (may lead to duplicate rows).
#' @param sf logical. Whether the output should be a spatial `"sf"` object 
#' using the Longitude and Latitude columns.
#' @importFrom readxl read_xlsx
#' @importFrom lubridate as_datetime
#' @importFrom dplyr rename mutate filter select full_join
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect starts_with
#' @importFrom sf st_as_sf
#' @importFrom magrittr "%>%"
#' @returns `list` containing the following objects:
#' \describe{
#' \item{`data`}{`"tbl_df"` object. `"sf"` if `sf == TRUE`.}
#' \item{`planar`}{Plane elements. Same row ids as in `data`.}
#' \item{`linear`}{Line elements. Same row ids as in `data`.}
#' }
#' @examples
#' \dontrun{
#' file = "E:/Lakehead/Field Work/StraboSpot_Output_11_01_2022_TS.xlsx"
#' dt <- read_strabo(file)
#' stereoplot(); 
#' stereo_greatcircle(dt$planar, col = "lightgrey") 
#' stereo_point(dt$linear)
#' }
read_strabo <- function(file, tag_cols = FALSE, sf = TRUE){
  data = readxl::read_xlsx(file, sheet = 1, skip = 2)
  colnames(data) <- make.names(colnames(data))
  
  data <- data %>% 
    dplyr::mutate(
      Date = lubridate::as_datetime(Date),
      Planar.Orientation.Dipdirection = (Planar.Orientation.Strike + 90) %% 360,
      Linear.Sense = ifelse(Planar.Orientation.Movement %in% c("left_lateral"), -1, NA),
      Linear.Sense = ifelse(Planar.Orientation.Movement  %in% c("right_lateral"), 1, Linear.Sense),
      Linear.Sense = ifelse(Planar.Orientation.Fault.Or.Sz.Type %in% c("sinistral", "reverse"), -1, Linear.Sense),
      Linear.Sense = ifelse(Planar.Orientation.Fault.Or.Sz.Type %in% c("dextral", "normal", "dextral_normal"), 1, Linear.Sense)
    )
  if(tag_cols){
    data <- data %>% 
    tidyr::pivot_longer(cols = tidyselect::starts_with("Tag"), names_to = "Tag", values_to = "u") %>%
    dplyr::filter(u == "X") %>%
    dplyr::select(-u) %>%
    dplyr::mutate(Tag = gsub("Tag:", "", Tag))
  }
  
  data.lines0 <- data %>%
    dplyr::select(tidyselect::starts_with("Linear")) %>%
    dplyr::filter(!is.na(Linear.Orientation.Trend))
  
  data.planes0 <- data %>%
    dplyr::select(tidyselect::starts_with("Planar"), ) %>%
    dplyr::filter(!is.na(Planar.Orientation.Dipdirection)) 
  
  
  data0 <- data %>% dplyr::select(
    !tidyselect::starts_with(c("Planar", "Linear"))
  ) 
  data0$Planar.Orientation.Unix.Timestamp = data$Planar.Orientation.Unix.Timestamp
  
  res <- dplyr::full_join(data.planes0, data.lines0,
              by = c("Planar.Orientation.Unix.Timestamp" = "Linear.Orientation.Unix.Timestamp"),
              multiple = "all"
    ) %>%
    dplyr::left_join(data0, by = "Planar.Orientation.Unix.Timestamp",
                     multiple = "all") %>%
    dplyr::mutate(Linear.Orientation.Unix.Timestamp = NA) %>%
    dplyr::select(colnames(data)) %>%
    dplyr::select(-Linear.Orientation.Unix.Timestamp)
  
  planes <- as.plane(cbind(res$Planar.Orientation.Dipdirection, res$Planar.Orientation.Dip))
  lines <- as.line(cbind(res$Linear.Orientation.Trend, res$Linear.Orientation.Plunge))
  
  if(sf){
    res = sf::st_as_sf(res, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
  } else {
    res = res
  }
  list(
    data = res,
    planar = planes,
    linear = lines
  )
}