#Step 1
#Functions for processing data inputs

#' Find the UTM zone for a given location
#'
#' This function returns the UTM zone for a given set of latitude/longitude.
#'
#' @usage zone_from_ll(lon, lat)
#'
#' @param lon A numeric vector of longitudes
#' @param lat A numeric vector of latitudes
#'
#' @details This function takes a point and calculates the correct UTM zone
#' in which that point falls. It was developed based on answers to a
#' Stack Overflow post (see below).
#'
#' @return Returns an integer vector of UTM zones the same length as the
#' inputs \code{lon} and \code{lat}
#'
#' @references Based on the answers to the Stack Overflow post found here: \cr
#' \url{https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude}
#'
#' @examples
#' #UTM grid for Ft. Lauderdale, FL, USA
#' zone_from_ll(lon = -80.15, lat = 26.133333)
#'
#' #UTM grid for Rome, Italy; Paris, France; London, UK
#' ll_df <- data.frame(city = c("Rome", "Paris", "London"),
#'                     lon = c(12.5, 2.3508, -0.1275),
#'                     lat = c(41.9, 48.8567, 51.507222))
#' zone_from_ll(lon = ll_df$lon, lat = ll_df$lat)
#'
#' @export
zone_from_ll <- function(lon, lat){
  #General form:
  zone <- (floor((lon + 180)/6) %% 60) + 1

  #Some rare exceptions to fix:
  zone2 <- case_when(
    lat > 55 & zone == 31 & lat < 64 & lon >2 ~ 32,
    lat > 71 & zone == 32 & lon < 9 ~ 31,
    lat > 71 & zone == 32 & lon > 8 ~ 33,
    lat > 71 & zone == 34 & lon < 21 ~ 33,
    lat > 71 & zone == 34 & lon > 20 ~ 35,
    lat > 71 & zone == 36 & lon < 33 ~ 35,
    lat > 71 & zone == 36 & lon > 32 ~ 37,
    TRUE ~ zone
  ) %>% as.integer()

  return(zone2)
}

#' Combine detections with stations
#'
#' This function combines detections with stations to place detections
#' in space.
#'
#' @usage proc_dets(det, sta, ...)
#'
#' @param det A \code{data.frame} of detections. Format should match the
#' object \strong{detections} from \link{acoustic}.
#' @param sta A \code{data.frame} of the stations where receivers were
#' located for a period of time. Format should match the object
#' \strong{stations} from \link{acoustic}.
#' @param ... Additional arguments (not currently implemented)
#'
#' @details Uses \code{\link[sqldf]{sqldf}} to combine \code{det} with
#' \code{sta} using \code{rec_id} where the date-time of the detection falls
#' within the date-time of the station deployment.
#'
#' It specifically excutes the following SQL \code{SELECT} query:
#'
#' \code{SELECT id, sta.rec_id, sta_id, dt, lon, lat} \cr
#' \code{FROM det LEFT JOIN sta}\cr
#' \code{ON det.rec_id = sta.rec_id}\cr
#' \code{AND det.dt > sta.dt_dep}\cr
#' \code{AND det.dt < sta.dt_ret}
#'
#' The left join keeps all rows from the detections
#' \code{data.frame}, inserting an \code{NA} if there is a
#' detection not within one of the time intervals for a given
#' receiver. The function subsequently removes any locations
#' with \code{NA} in the x- or y-coordinate with a warning.
#'
#' @return Returns a \code{data.frame} of georeferenced detections.
#'
#' @seealso \link{acoustic} for details on the formatting of the
#' \code{data.frame}s \code{det} and \code{sta}
#'
#' @examples
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' @export
proc_dets <- function(det, sta){
  #Check that the dt_* columns are formatted
  if((!lubridate::is.POSIXct(sta$dt_dep) | !lubridate::is.POSIXct(sta$dt_ret))){
    stop("The dt columns of 'sta' should be formatted as POSIXct.")
  }
  if(!lubridate::is.POSIXct(det$dt)){
    stop("The dt column of 'det' should be formatted as POSIXct.")
  }
  #Use SQL join to join 'det' and 'sta'
  res <- sqldf::sqldf("SELECT id, sta.rec_id, sta_id,
                        dt, x, y
                      FROM det LEFT JOIN sta
                        ON det.rec_id = sta.rec_id
                          AND det.dt > sta.dt_dep
                          AND det.dt < sta.dt_ret")
  #Remove any detections that were not assigned to a location
  rem <- which(is.na(res$x) | is.na(res$y))
  if (length(rem) > 0) {
    #Report to the user.
    warning(paste0("\n", length(rem), " out of ", nrow(res),
                " detections were not successfully associated with a station.",
                "\n\nThese detections have been removed."))
    #Drop those rows
    res <- res[-rem,]
  }

  return(res)
}

#' Plots processed station history
#'
#' This function takes a processed detection history returned by
#' \code{\link{proc_dets}()} and plots the station history.
#'
#' @usage plot_sta_history(proc_det, set.par = TRUE, ...)
#'
#' @param proc_det A \code{data.frame} of georeferenced detections as returned
#' by the function \code{\link{proc_dets}()}.
#' @param set.par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? (This should be \code{FALSE} only if the
#' user wishes to manually set the graphical parameters).
#' @param ... Additional arguments (not currently implemented)
#'
#' @details Details here
#'
#' @return Plots a figure of station activity over time.
#'
#' @seealso \code{\link{proc_dets}} for details on the formatting of the
#' \code{data.frame} \code{proc_det}
#'
#' @examples
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Plot station history
#' plot_sta_history(proc.det)
#'
#' @export
plot_sta_history <- function(proc_det, set.par=TRUE, ...) {
  #Store original parameters
  orig.par <- par(no.readonly = TRUE)
  #Change parameters
  if(set.par){
    par(mar = c(2.5, 6.5, 2.0, 1.0),
        cex = 0.8)
  }
  #Factor the station id
  proc_det$sta_id <- factor(proc_det$sta_id)
  levs <- levels(proc_det$sta_id)
  #Plot
  {plot(x = proc_det$dt, y = proc_det$sta_id, pch = 16,
       axes = FALSE, xlab=NA, ylab=NA, main="Detections by Station")
  axis.POSIXct(1, proc_det$dt)
  axis(2, labels=levs, at=1:length(levs), las=2)
  box()}
  #Revert to original parameters
  if(set.par){
    on.exit(par(orig.par), add = TRUE)
  }
}

#' Plots a map of station locations
#'
#' This function takes a processed detection history returned by
#' \code{\link{proc_dets}()} and plots a map showing the relative number
#' of detections per station.
#'
#' @usage map_dets(proc_det, sta.crs = 4326, base.layers = NULL,
#'              base.borders = NULL, base.cols = NULL,
#'              sta.col = "black", sta.bg = "red",
#'              leg.pos = "bottomleft", set.par = TRUE, ...)
#'
#' @param proc_det A \code{data.frame} of georeferenced detections as returned
#' by the function \code{\link{proc_dets}()}.
#' @param sta.crs A coordinate reference system to use for the station
#' locations. Used by \code{\link[sf]{st_sf}} as the \code{crs} argument.
#' Defaults to integer 4326, indicating EPSG code 4326 for WGS 84 lat/long.
#' @param base.layers A \code{list} of georeferenced vectors to be plotted
#' as the basemap under the stations. Designed to be \code{\link[sf]{sfc}}
#' \code{LINESTRING} or \code{POLYGON} objects.
#' @param base.borders A \code{list} of colors to used to outline the
#' polygons listed in \code{base.layers}. Should be the same length as
#' \code{base.layers}.
#' @param base.cols A \code{list} of colors to used to fill the polygons
#' listed in \code{base.layers}.  Should be the same length as \code{
#' base.layers}.
#' @param leg.pos Indicates where the legend should be. Can be any \code{x}
#' that the base \code{graphics} function \code{\link[graphics]{legend}}
#' will accept. Defaults to the text string \code{"bottomleft"}.
#' @param set.par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? This should be \code{FALSE} if the
#' user wishes to manually set the graphical parameters, \emph{e.g.}, so
#' that they can still add to the plot manually after it is generated.
#' @param ... Additional arguments (not currently implemented)
#'
#' @details Details here
#'
#'
#' @return Plots a figure of station activity over time.
#'
#' @seealso \code{\link{proc_dets}} for details on the formatting of the
#' \code{data.frame} \code{proc_det}
#'
#' @examples
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Map detections with all defaults
#' map_dets(proc.det)
#'
#' #Map detections with custom settings
#' map_dets(proc_det = proc.det, sta.crs = 4326,
#'     base.layers = list(acoustic$study_area, acoustic$land),
#'     base.cols = list("gray30", "wheat"), base.borders = list(NA, "seagreen"),
#'     xlim=c(-64.628, -64.612), ylim=c(17.770, 17.795),
#'     leg.pos = "topleft", sta.col = "blue", sta.bg = "orange")
#'
#' @export
map_dets <- function(proc_det, sta.crs = 4326, base.layers = NULL,
                     base.borders = NULL, base.cols = NULL,
                     sta.col = "black", sta.bg = "red",
                     leg.pos = "bottomleft",
                     xlim = NULL, ylim = NULL, set.par = TRUE, ...) {
  #Store original parameters
  orig.par <- par(no.readonly = TRUE)
  #Change parameters
  if(set.par){
    par(mar = c(0.1, 0.1, 0.1, 0.1),
        cex = 0.8)
  }

  #Aggregate station detections
  stas <- proc_det %>%
    group_by(sta_id) %>%
    summarize(n_det = n()) %>%
    left_join(unique(proc_det[,c("sta_id", "x", "y")]), by = "sta_id") %>%
    mutate(sz = cut(n_det, breaks=pretty(n_det), labels=paste0("<", pretty(n_det)[-1])))

  #Create an sf object
  sta.sp <- sf::st_geometry(
    sf::st_as_sf(stas, coords = c("x", "y"),
                 dim = "XY", crs = sta.crs))

  #Find xlim and ylim by checking extent of all layers
  #Combine all layers in a list
  if(is.null(base.layers)){
    all.layers <- list()
  } else {
    all.layers <- base.layers
  }
  all.layers[[length(base.layers)+1]] <- sta.sp

  if (is.null(xlim)){  #If xlim is not specified, pull from layers
    #xmin
    xmin <- min(unlist(lapply(all.layers, function(x) st_bbox(x)["xmin"])))
    #xmax
    xmax <- max(unlist(lapply(all.layers, function(x) st_bbox(x)["xmax"])))
    #Combine
    xlim <- c(xmin, xmax)
  }

  if (is.null(ylim)){  #If ylim is not specified, pull from layers
    #ymin
    ymin <- min(unlist(lapply(all.layers, function(x) st_bbox(x)["ymin"])))
    #ymax
    ymax <- max(unlist(lapply(all.layers, function(x) st_bbox(x)["ymax"])))
    #Combine
    ylim <- c(ymin, ymax)
  }

  #Plot blank basemap
  #Create sf point object
  dummy.point <- sf::st_point(x = c(mean(xlim), mean(ylim)), dim="XY")
  plot(dummy.point, xlim = xlim, ylim = ylim, pch=NA)

  #Add to basemap first if base.layers is not NULL
  if(!is.null(base.layers)){
    #First check that the color lists are the right lengths
    ##Borders
    if (!is.null(base.borders)){
      #Make sure it is same length as base.layers
      if(length(base.layers) != length(base.borders)){
        stop("If you specify 'base.borders', the list must be the
               same length as 'base.layers'. Note that you can pad this
               list with NAs where necessary.")
      }
    } else {
      base.borders <- rep(NA, length(base.layers))
    }
    ##Colors
    if(!is.null(base.cols)){
      #Make sure it has the same length as base.layers
      if(length(base.layers) != length(base.cols)){
        stop("If you specify 'base.cols', the list must be the
               same length as 'base.layers'. Note that you can pad this
               list with NAs where necessary.")
      }
    } else {
      base.cols <- rep(NA, length(base.layers))
    }
    #Now that we have the colors, plot each layer
    #Loop through each layer
    for (i in 1:length(base.layers)){

      #Should be either "sfc_LINESTRING" or "sfc_POLYGON"
      if(!any(class(base.layers[[i]]) %in%
              c("sfc_LINESTRING", "sfc_POLYGON"))){

        stop(paste("Element", i, "of list 'base.layers' is neither an
                    'sfc_LINESTRING' nor 'sfc_POLYGON'."))

      } else { #As long as the classes are right

        #If it's a LINESTRING
        if (any(class(base.layers[[i]])=="sfc_LINESTRING")){
          plot(base.layers[[i]], add=TRUE,
               col=ifelse(is.na(base.cols[[i]]), 1, base.cols[[i]]))
        }

        #If it's a POLYGON
        if (any(class(base.layers[[i]])=="sfc_POLYGON")){
          plot(base.layers[[i]], add=TRUE,
               col=base.cols[[i]],
               border=ifelse(is.na(base.borders[[i]]), 1, base.borders[[i]]))
        }
      }
    }
  } #End of basemap plotting

  #Plot stations
  plot(sta.sp, add = TRUE,
       pch = 21, col = sta.col, bg = sta.bg,
       cex = as.numeric(stas$sz)/2)
  #Plot legend
  legend(x = leg.pos, legend = levels(stas$sz),
         pch = 21, col = sta.col, pt.bg = sta.bg,
         pt.cex = (1:length(levels(stas$sz)))/2,
         title = expression(bold("Detections")))


  #Revert to original parameters
  if(set.par){
    on.exit(par(orig.par), add = TRUE)
  }
  #Return the aggregated data
  return(stas)
}
