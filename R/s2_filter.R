#Step 2
#Functions for filtering possibly erroneous detections

#Speed filter----


spd_filter <- function(dets, max_speed, rec_range = 500){
  #Check class
  if(!is.dets(dets)){
    stop("Object 'dets' must be of class 'dets'. See ?proc_dets for details.")
  }

  #Need to convert from lat/long to UTM
  if(grepl("longlat", sf::st_crs(paths$geometry)$proj4string)){
    warning(paste("Function must convert from a geographic coordinate",
                  "system to a projected coordinate system. This is safest",
                  "if the original CRS was EPSG:4326, i.e., st_crs(4326)."))
    #Store original CRS
    orig_crs <- sf::st_crs(paths$geometry)
    #Get UTM zone
    mean_coords <- paths %>%
      summarize(x = mean(x1, na.rm = TRUE),
                y = mean(y1, na.rm = TRUE))
    z <- zone_from_ll(lon = mean_coords[["x"]],
                      lat = mean_coords[["y"]])
    utm_crs <- as.integer(paste0(326, z))
    #Transform
    paths$geometry <- sf::st_transform(paths$geometry, crs = utm_crs)
  }

  #Store projected CRS
  proj_crs <- sf::st_crs(paths$geometry)


}
