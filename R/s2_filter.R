#Step 2
#Functions for filtering possibly erroneous detections

#Speed filter----

#' Filters out detections based on maximum speed over distance
#'
#' Filters detections that would imply an unreasonably high swim speed
#' over a sufficiently large distance.
#'
#' @param proc_det A \code{data.frame} of class \code{dets} as returned
#' by the function \code{\link{proc_dets}()}.
#' @param max_speed An object of class \code{\link[units]{units}} setting
#' the maximum speed with which an animal can move between two receivers that
#' are at least \code{min_dist} apart.
#' @param min_dist An object of class \code{\link[units]{units}} setting the
#' minimum distance two receivers must be apart to consider the speed filter.
#'
#' @details False-positive detections in acoustic data can occur for reasons
#' related to signal collisions. Type A false detections are detections of
#' IDs that are not present in the study system, and here we assume the user
#' has already removed these detections prior to using this package. Type B
#' false detections occur when a signal is erroneously decoded as the signal
#' of a different tag that is in the study system -- these detections are
#' much harder to distinguish from true detections and are the focus of this
#' function.
#'
#' Determining the maximum speed an animal can move is not as
#' straightforward as it first seems. Animals can move at very high speeds
#' for very short bursts, but most often cannot sustain those high speeds
#' over longer distances. Also, with the considerable position error that
#' can result from assuming a detected animal is at the location of the
#' acoustic station, an animal sitting in the middle of two receivers and
#' being detected simultaneously by both can seem like an animal moving with
#' impossible speed.
#'
#' @seealso See Simpfendorfer \emph{et al.}
#' (\href{https://doi.org/10.1186/s40317-015-0094-z}{2015})
#' for details about false detections in acoustic telemetry.
#'
#' @examples
#'
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Apply the speed filter with default max_speed and min_dist
#' det.filt <- spd_filter(proc.det)
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

#Single detection filter----
#' Filters out singleton detections
#'
#' Filters out detections if only 1 detection occurred for an ID in
#' a user-defined time horizon.
#'
#' @usage singleton_filter(proc_det, time_horizon = "1 day")
#'
#' @param proc_det A \code{data.frame} of class \code{dets} as returned
#' by the function \code{\link{proc_dets}()}.
#' @param time_horizon The time during which a single detection will be
#' filtered. Can be any value that can be accepted by the \code{by} argument
#' of \code{\link{seq.POSIXt}()}.
#'
#' @details False-positive detections in acoustic data can occur for reasons
#' related to signal collisions. Type A false detections are detections of
#' IDs that are not present in the study system, and here we assume the user
#' has already removed these detections prior to using this package. Type B
#' false detections occur when a signal is erroneously decoded as the signal
#' of a different tag that is in the study system -- these detections are
#' much harder to distinguish from true detections and are the focus of this
#' function.
#'
#' This function is based on the assumption that it is possible for signal
#' collisions to result in a single false detection over a sufficiently long
#' time horizon, but that the probability of two false detections of the same
#' ID over that same time horizon is very low. Thus, the choice of time horizon
#' is very critical here. The period of time should be long enough that your
#' animal would plausibly use a receiver many times, but not so long that the
#' probability of two false detections of the same ID becomes unacceptably
#' high.
#'
#' @seealso See Simpfendorfer \emph{et al.}
#' (\href{https://doi.org/10.1186/s40317-015-0094-z}{2015})
#' for details about false detections in acoustic telemetry.
#'
#' @examples
#'
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Apply the singleton filter with default time_horizon
#' det.filt <- singleton_filter(proc.det)
#'
#' #Apply the singleton filter with 1 month time horizon
#' ##Note: The user should see ?lubridate::add_with_rollback
#'  #to understand the behavior when adding months.
#' det.filt2 <- singleton_filter(proc.det, time_horizon = "1 month")
#'
#' @export
singleton_filter <- function(proc_det, time_horizon = "1 day"){

  #Check class
  if(!is.dets(proc_det)){
    stop("Object 'proc_det' must be of class 'proc_det'. See ?proc_proc_det for details.")
  }

  #Unique ids
  uids <- unique(proc_det$id)

  #Store CRS
  crs <- sf::st_crs(proc_det)

  #Initialize data.frame of flagged detections
  flagged_dets <- data.frame()

  #Report to user
  cat(paste0("Filtering singletons from ", length(uids), " IDs...\n"))

  #Loop through each ID
  for(i in 1:length(uids)){
    #Report to user
    cat(paste0("  ... ID: ", uids[i], "\n"))
    #Subset detections
    dets_sub <- proc_det %>%
      sf::st_drop_geometry() %>%
      filter(id == uids[i])
    #Get endpoints
    endpoints <- dets_sub %>%
      arrange(dt) %>%
      summarize(start_dt = min(dt, na.rm = TRUE),
                end_dt = max(dt, na.rm = TRUE)) %>%
      mutate(start_dt = as.POSIXct(trunc(start_dt, "hours")),
             end_dt = as.POSIXct(trunc(end_dt, "hours") +
                                   lubridate::hours(1))) %>%
      #Need add_with_rollback in the case the user adds a period with months
      mutate(end_dt = lubridate::add_with_rollback(end_dt, lubridate::period(time_horizon)))
    #Create intervals
    ints <- ADePTR:::s_e_to_interval(endpoints, time_horizon = time_horizon)
    #Create interval data.frame
    int_df <- data.frame(int_id = 1:length(ints),
                         ints = ints)
    #Assign interval id to each location
    dets_sub <- dets_sub %>%
      rowwise() %>%
      mutate(int_id = int_df[which(lubridate::`%within%`(dt, int_df$ints)), "int_id"]) %>%
      ungroup()
    #Identify those intervals which only have 1 detection at a station
    singletons <- dets_sub %>%
      group_by(int_id, sta_id) %>%
      tally() %>%
      mutate(flag = case_when(
        n == 1 ~ TRUE,
        TRUE ~ FALSE))
    #Left join to dets_sub
    dets_flagged <- dets_sub %>%
      left_join(singletons, by = c("int_id", "sta_id"))
    #Add to flagged_dets
    flagged_dets <- bind_rows(flagged_dets, dets_flagged)
  }

  #Separate into flagged and good
  flagged <- flagged_dets %>%
    filter(flag)
  good <- flagged_dets %>%
    filter(!flag) %>%
    dplyr::select(id, rec_id, sta_id, dt, x, y)

  #Report to the user
  cat(paste0("  ... Done. \n", nrow(flagged), " detections were removed.\n"))

  #Restore geometry
  good <- good %>%
    sf::st_as_sf(agr = "identity", coords = c("x", "y"), crs = crs, remove = FALSE)

  #Should still be same class as proc_det
  class(good) <- class(proc_det)

  #Return
  return(good)
}


s_e_to_interval <- function(x, time_horizon){
  #Starts
  starts <- seq(from = x[["start_dt"]],
                to = x[["end_dt"]],
                by = time_horizon)
  #Ends
  ends <- lead(starts) - lubridate::seconds(1)
  #Intervals
  ints <- lubridate::interval(start = starts, end = ends)
  #Drop the last one
  ints <- ints[-length(ints)]
  #Return
  return(ints)
}
