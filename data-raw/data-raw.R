#This script processes raw data into RData formats that can be included
  #in the package.
library(dplyr)
library(sf)
library(lubridate)

#Stations----
stations <- read.csv("data-raw/dwnld_metadata_all_updtd01022018.csv", stringsAsFactors = FALSE) %>%
  select(sta_id, rec_id, dt_dep=dep_dt, dt_ret=ret_dt, x=sta_lon, y=sta_lat) %>%
  filter(sta_id %in% sprintf("BUIS_%02d", seq(1, 52, by=1))) %>%
  mutate(dt_dep=mdy_hm(dt_dep), dt_ret=mdy_hm(dt_ret))

sta_loc <- stations %>%
  select(sta_id, x, y) %>%
  st_as_sf(coords=c("x", "y"), crs=4326)

#Buck Island----
bi <- sf::st_geometry(
  sf::st_transform(
    sf::st_read("data-raw/Buck_Island_EPSG26920.shp"), 4326))

#BIRNM----
birnm <- sf::st_geometry(
  sf::st_transform(
    sf::st_read("data-raw/birnm_2011boundary_prjctd.shp"), 4326))


#Turtle detections----
hawk <- read.csv("data-raw/hawk_data_final_04212016.csv", stringsAsFactors = FALSE) %>%
  select(t_id, t_name, rec_id, dt=dt_lt, sta_id) %>%
  mutate(dt = ymd_hms(dt), ym = format(dt, "%Y-%m")) %>%
  filter(ym %in% paste0("2014-", c("05", "06", "07", "08"))) %>%
  #filter(t_name %in% c(44455, 32174)) %>%
  filter(t_name %in% c(44445, 32174, 44455, 44452)) %>%
  filter(sta_id %in% stations$sta_id) %>%
  select(id=t_id, dt, rec_id)

#Make object to save----
acoustic <- list(stations=stations,
                 detections=hawk,
                 land=bi,
                 study_area=birnm)

#Save to package----
usethis::use_data(acoustic, overwrite = TRUE)
