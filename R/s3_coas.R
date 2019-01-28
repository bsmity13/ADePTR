#Step 3
#Calculate Centers of Activity (COAs) from detections

#' Calculate Harmonic Mean
#'
#' This function is used internally by \code{\link{coa_locs}()} to
#' calculate the harmonic mean of a set of detections.
#'
#' @usage harmonic_mean(x)
#'
#' @param x A vector of numbers
#'
#' @details Calculates the harmonic mean
#' (\href{https://en.wikipedia.org/wiki/Harmonic_mean}{Wikipedia}), which
#' is the reciprocal of the arithmetic mean of the reciprocals, i.e.:
#'
#' \code{1/(mean(1/x))}
#'
#' @return Returns a vector of length 1
#'
#' @seealso Used internally by \code{\link{coa_locs}}
harmonic_mean <- function(x){
  return(1/(mean(1/x)))
}

#' Calculate Centers of Activity
#'
#' This function calculates Centers of Activity (COAs) by averaging
#' the locations of detections in a user-defined time period.
#'
#' @usage coa_locs(proc_det, Delta_t = "1 hour", mean_type = c("arithmetic",
#'  "harmonic"), ...)
#'
#' @param proc_det A \code{data.frame} of georeferenced detections as returned
#' by the function \code{\link{proc_dets}()}.
#' @param Delta_t The desired time interval for the COAs (\eqn{\Delta}t).
#' @param mean_type The type of mean used to calculate the mean positions.
#' Should be either \code{"arithmetic"} or \code{"harmonic"}.
#' @param ... Additional arguments (not currently implemented)
#'
#' @details This function implements the two mean position algorithms
#' suggested by Simpfendorfer et al. (2002) for calculating centers of
#' activity (COAs).
#'
#' @return Returns a \code{data.frame} of georeferenced centers of activity.
#'
#' @seealso \code{\link{proc_dets}} for details on the formatting of the
#' \code{data.frame} \code{proc_det}
#'
#' See Simpfendorfer \emph{et al.} (\href{https://doi.org/10.1139/f01-191}{2002}) for details on using COAs to study
#' marine animal movements.
#'
#' @references
#'
#' Simpfendorfer, C.A., M.R. Heupel, and R.E. Hueter. 2002. Estimation of
#' short-term centers of activity from an array of omnidirectional hydrophones
#' and its use in studying animal movements. \emph{Canadian Journal of
#' Fisheries and Aquatic Sciences} 59(1): 23-32.
#'
#' @export
coa_locs <- function(proc_det, Delta_t = "1 hour", mean_type = c("arithmetic", "harmonic"), ...){

  #Check that mean_type is one of the possible values
  if(!(mean_type[1] %in% c("arithmetic", "harmonic"))){
    stop("'mean_type' must be either 'arithmetic' or 'harmonic'")
  }

  #Find the start time for all COAs and round down to the hour
  s.time <- lubridate::ymd_hms(
    format(min(proc_det$dt), "%Y-%m-%d %H:00:00"), tz="")
  #Find the end time for all COAs and round up to the next hour
  e.time <- lubridate::ymd_hms(
    format(max(proc_det$dt), "%Y-%m-%d %H:00:00"), tz="") +
    lubridate::hours(1)

  #Create vector of interval start times using Delta_t
  int.starts <- seq(s.time, e.time, by=Delta_t)
  #Creat vector of interval end times
  int.ends <- int.starts - lubridate::seconds(1)
  #Drop the last start and the first end
  int.starts <- int.starts[-length(int.starts)]
  int.ends <- int.ends[-1]
  #Combine in data.frame of intervals
  ints <- data.frame(
    s = int.starts,
    e = int.ends,
    int = as.character(lubridate::interval(int.starts, int.ends)),
    int_num = 1:length(int.starts)
  )

  #Join proc_det with ints using SQL
  det_ints <- sqldf::sqldf("SELECT id, dt, x, y, int, int_num
               FROM proc_det
                  LEFT JOIN
                  ints ON proc_det.dt >= ints.s
                      AND proc_det.dt <= ints.e")

  #Calculate mean location
  if(mean_type[1] == "arithmetic"){
    coas <- det_ints %>%
      group_by(id, int_num) %>%
      summarize(x = mean(x), y = mean(y), dt = mean(dt), locs = n()) %>%
      ungroup() %>%
      select(id, dt, x, y, int_num, locs) %>%
      arrange(id, dt) %>%
      as.data.frame
  }

  if(mean_type[1] == "harmonic"){
    coas <- det_ints %>%
      group_by(id, int_num) %>%
      summarize(x = harmonic_mean(x), y = harmonic_mean(y), dt = mean(dt), locs = n()) %>%
      select(id, dt, x, y, int_num, locs) %>%
      arrange(id, dt) %>%
      as.data.frame
  }
return(coas)
}

#' Map Centers of Activity
#'
#' Plots a map using \code{\link{map_dets}()} and adds the
#' COAs returned by \code{\link{coa_locs}()} to the map.
#'
#' @usage
#' map_coas(proc_det, coas, coa.crs = 4326, coa.palette = viridis::viridis,
#' coa.leg.pos = "bottomleft", det.leg.pos = NA, set.par = TRUE, ...)
#'
#' @param proc_det The \code{data.frame} of processed detections used by
#' \code{\link{coa_locs}()} to create the COAs.
#' @param coas The \code{data.frame} of centers of activiy (COAs) created
#' by \code{\link{coa_locs}()}.
#' @param coa.crs A coordinate reference system to use for the COA
#' locations. Used by \code{\link[sf]{st_sf}} as the \code{crs} argument.
#' Defaults to integer 4326, indicating EPSG code 4326 for WGS 84 lat/long.
#' @param coa.palette A function name that will generate a palette for each
#' of the individuals tracked. Can be any function that will accept \code{n}
#' as an argument.
#' @param coa.leg.pos Indicates where the COA legend should be (showing
#' the symbol for each individual animal). Can be any \code{x}
#' that the base \code{graphics} function \code{\link[graphics]{legend}}
#' will accept. Defaults to the text string \code{"bottomleft"}.
#' @param det.leg.pos Indicates where the default legend for the \emph{
#' station detections} should be. Passed to \code{\link{map_dets}()} to
#' control that legend placement. Can be any \code{x} that the base
#' \code{graphics} function \code{\link[graphics]{legend}} will accept.
#' Defaults to \code{NA} which prevents that legend from being drawn. Can
#' be changed to a location \code{!= coa.leg.pos}.
#' @param set.par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? This should be \code{FALSE} if the
#' user wishes to manually set the graphical parameters, \emph{e.g.}, so
#' that they can still add to the plot manually after it is generated.
#' @param ... Additional arguments to be passed to \code{\link{map_dets}},
#' \emph{e.g.}, \code{base.layers}, \code{xlim}, \code{ylim}, \emph{etc}.
#'
#' @details Calls \code{\link{map_dets}()} to draw the basemap and then adds
#' the COAs on top.
#'
#' @seealso \code{\link{map_dets}}
#'
#' @export
map_coas <- function(proc_det, coas, coa.crs = 4326, coa.palette = viridis::viridis,
                     coa.leg.pos = "bottomleft", det.leg.pos = NA,
                     set.par = TRUE, ...){

  #Store original parameters
  orig.par <- par(no.readonly = TRUE)
  #Change parameters
  if(set.par){
    par(mar = c(0.1, 0.1, 0.1, 0.1),
        cex = 0.8)
  }

  #Factor id
  coas$id <- factor(coas$id)

  #Create an sfc_POINT object
  coas.sp <- sf::st_as_sf(coas, coords = c("x", "y"),
                          dim = "XY", crs = coa.crs) %>%
    sf::st_geometry()

  #Create an sfc_MULTILINESTRING object
  coas.lines <- coas %>%
    sf::st_as_sf(coords = c("x", "y"), dim = "XY", crs = coa.crs) %>%
    group_by(id) %>%
    summarize() %>%
    sf::st_cast("MULTILINESTRING") %>%
    sf::st_geometry()

  #Pass proc_det and ... to map_dets
  map_dets(proc_det = proc_det, leg.pos = det.leg.pos, set.par = FALSE, ...)

  #Setup palette for individuals
  orig.pal <- palette() #Get previous palette
  n.ind <- length(unique(coas$id))
  palette(do.call(coa.palette, args = list(n = n.ind))) #Change palette

  #Add lines to map
  plot(coas.lines, col = 1:length(levels(coas$id)), add=TRUE)

  #Add COAs to map
  plot(coas.sp, pch = 22, cex = 0.8,
       col = "black", bg = as.numeric(coas$id),
       add = TRUE)

  #Plot legend
  legend(x = coa.leg.pos, legend = levels(coas$id),
         pch = 22, pt.cex = 0.8, col = "black",
         pt.bg = 1:length(levels(coas$id)),
         title = expression(bold("Individuals")))

  #Reset palette
  palette(orig.pal)


  #Revert to original parameters
  if(set.par){
    on.exit(par(orig.par), add = TRUE)
  }

}
