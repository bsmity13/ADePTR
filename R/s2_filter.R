#Step 2
#Functions for filtering possibly erroneous detections

#Speed filter----

#' Filters out detections based on maximum speed over distance
#'
#' DOCUMENT ME, PLEASE
#'
#' @export
spd_filter <- function(proc_det,
                       max_speed = units::set_units(10, "m/s"),
                       min_dist = units::set_units(1, "km")){
  #Check class
  if(!is.dets(proc_det)){
    stop("Object 'proc_det' must be of class 'proc_det'. See ?proc_proc_det for details.")
  }

  #Need to convert from lat/long to UTM
  if(grepl("longlat", sf::st_crs(proc_det)$proj4string)){
    warning(paste("Function must convert from a geographic coordinate",
                  "system to a projected coordinate system. This is safest",
                  "if the original CRS was EPSG:4326, i.e., st_crs(4326)."))
    #Store original CRS
    orig_crs <- sf::st_crs(proc_det)
    #Get UTM zone
    mean_coords <- proc_det %>%
      summarize(x = mean(x, na.rm = TRUE),
                y = mean(y, na.rm = TRUE))
    z <- zone_from_ll(lon = mean_coords[["x"]],
                      lat = mean_coords[["y"]])
    utm_crs <- as.integer(paste0(326, z))
    #Transform
    proc_det <- sf::st_transform(proc_det, crs = utm_crs)
    proc_det$x <- sf::st_coordinates(proc_det)[,1]
    proc_det$y <- sf::st_coordinates(proc_det)[,2]
  }

  #Store projected CRS
  proj_crs <- sf::st_crs(proc_det)

  #Calculate speed
  spd <- proc_det %>%
    arrange(id, dt) %>%
    mutate(x1 = x,
           y1 = y,
           x2 = lead(x, n = 1L),
           y2 = lead(y, n = 1L),
           dt2 = lead(dt, n = 1L)) %>%
    mutate(meters = units::set_units(
      round(sqrt((x1 - x2)^2 + (y1 - y2)^2), 1), "m"),
           secs = units::set_units(
             as.numeric(
             difftime(dt2, dt, units = "secs")), "s")) %>%
    mutate(speed = case_when(
      meters < units::set_units(1, "m") ~ units::set_units(0, "m/s"),
      TRUE ~ round(meters/secs, 2)
    ))

  #Filter if speed parameters are exceeded
  spd_filt <- spd %>%
    mutate(flag = case_when(
      meters >= min_dist & speed > max_speed ~ TRUE,
      TRUE ~ FALSE
    ))

  #Report to the user and remove
  cat(paste0(sum(spd_filt$flag), " detections were removed.\n"))
  spd_filt <- spd_filt %>%
    filter(!flag)

  #Convert to orginal CRS
  if(sf::st_crs(spd_filt) != orig_crs){
    spd_filt <- sf::st_transform(spd_filt, orig_crs) %>%
      mutate(x = sf::st_coordinates(geometry)[,1],
             y = sf::st_coordinates(geometry)[,2])
  }

  #Keep necessary columns
  spd_filt <- spd_filt %>%
    dplyr::select(id, rec_id, sta_id, dt, x, y, geometry)

  #Should still be same class as proc_det
  class(spd_filt) <- class(proc_det)

  return(spd_filt)
}
