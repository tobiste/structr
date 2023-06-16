#' Import orientation data from *StraboSpot* field book
#' 
#' Import the XLS export of field book data or StraboMobile data from 
#' \url{strabospot.org/my_data} into R's work space
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
#' \item{`planar`}{Plane elements. Same row IDs as in `data`.}
#' \item{`linear`}{Line elements. Same row IDs as in `data`.}
#' }
#' @name strabo
#' @examples
#' \dontrun{
#' file = "E:/Lakehead/Field Work/StraboSpot_Output_11_01_2022_TS.xlsx"
#' dt <- read_strabo(file)
#' stereoplot(); 
#' stereo_greatcircle(dt$planar, col = "lightgrey") 
#' stereo_point(dt$linear)
#' 
#' read_strabo_mobile("C:/Users/tobis/Downloads/StraboSpot_Search_06_16_2023.txt")
#' }
NULL

#' @rdname strabo
#' @export
read_strabo <- function(file, tag_cols = FALSE, sf = TRUE){
  Date <- Planar.Orientation.Dipdirection <- Planar.Orientation.Strike <- 
    Linear.Sense <- Planar.Orientation.Movement <- temp <- Tag <- Linear.Orientation.Trend <- 
    Planar.Orientation.Fault.Or.Sz.Type <- Planar.Orientation.Unix.Timestamp <- Linear.Orientation.Unix.Timestamp <- 
    Planar.Orientation.Dip <- Linear.Orientation.Plunge <- Longitude <- Latitude <- NULL
  
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
    tidyr::pivot_longer(cols = tidyselect::starts_with("Tag"), names_to = "Tag", values_to = "temp") %>%
    dplyr::filter(temp == "X") %>%
    dplyr::select(-temp) %>%
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

#' @rdname strabo
#' @export
read_strabo_mobile <- function(file, sf = TRUE){
  data0 <- read.table(x, header = TRUE, sep = "\t", colClasses = c("integer", rep("character", 3), rep("numeric", 11), "character")) 
  
  lines0 <- dplyr::filter(data0, Type == "L")
  lines <- as.line(cbind(lines0$Trd.Strk, lines0$Plg.Dip))
  lines.meta <- lines0 %>% dplyr::select(-c(No., Trd.Strk, Plg.Dip))
  rownames(lines) <- rownames(lines.meta) <- lines0$No.
  
  planes0 <-  dplyr::filter(data0, Type == "P") %>%
    mutate(Dipdir = (Trd.Strk + 90) %% 360)
  planes <- as.plane(cbind(planes0$Dipdir, planes0$Plg.Dip))
  planes.meta <- planes0 %>% dplyr::select(-c(No., Dipdir, Trd.Strk, Plg.Dip))
  rownames(planes) <- rownames(planes.meta) <- planes0$No.
  
  
  if(sf){
    lines.meta = sf::st_as_sf(lines.meta, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
    planes.meta = sf::st_as_sf(planes.meta, coords = c("Longitude", "Latitude"), crs = "WGS84", remove = FALSE, na.fail = FALSE)
  } 
  
  list(
    linear = lines,
    linear_data = lines.meta,
    planar = planes,
    planar_data = planes.meta
  )
}