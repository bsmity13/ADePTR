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
#' @usage proc_dets(det, sta, crs = 4326)
#'
#' @param det A \code{data.frame} of detections. Format should match the
#' object \code{detections} from \link{acoustic} (\code{id}, \code{dt},
#' \code{rec_id}).
#' @param sta A \code{data.frame} of the stations where receivers were
#' located for a period of time. Format should match the object
#' \code{stations} from \link{acoustic} (\code{sta_id}, \code{rec_id},
#' \code{dt_dep}, \code{dt_ret}, \code{x}, \code{y}).
#' @param crs Coordinate Reference System to use for the detections. Passed
#' to \code{\link[sf:st_crs]{sf::st_crs}()} to set CRS for sf object.
#' Defaults to \code{4326}, longitude/latitude on the WGS84 spheroid.
#'
#' @details Uses \code{\link[sqldf]{sqldf}} to combine \code{det} with
#' \code{sta} using \code{rec_id} where the date-time of the detection falls
#' within the date-time of the station deployment.
#'
#' It specifically excutes the following SQL \code{SELECT} query:
#'
#' \preformatted{
#'  SELECT id, sta.rec_id, sta_id, dt, lon, lat
#'  FROM det LEFT JOIN sta
#'  ON det.rec_id = sta.rec_id
#'      AND det.dt > sta.dt_dep
#'      AND det.dt < sta.dt_ret}
#'
#' The left join keeps all rows from the detections
#' \code{data.frame}, inserting an \code{NA} if there is a
#' detection not within one of the time intervals for a given
#' receiver. The function subsequently removes any locations
#' with \code{NA} in the x- or y-coordinate with a warning.
#'
#' @return Returns a \code{data.frame} of class \code{dets} of
#' georeferenced detections. Object is also of class \code{sf} and
#' contains column \code{geometry}.
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
proc_dets <- function(det, sta, crs = 4326){
  #Check that the dt_* columns are formatted
  if((!lubridate::is.POSIXct(sta$dt_dep) | !lubridate::is.POSIXct(sta$dt_ret))){
    stop("The dt columns of 'sta' should be formatted as POSIXct.")
  }
  if(!lubridate::is.POSIXct(det$dt)){
    stop("The dt column of 'det' should be formatted as POSIXct.")
  }

  cat("\nJoining data... ")
  #Use SQL join to join 'det' and 'sta'
  res <- sqldf::sqldf("SELECT id, sta.rec_id, sta_id,
                        dt, x, y
                      FROM det LEFT JOIN sta
                        ON det.rec_id = sta.rec_id
                          AND det.dt > sta.dt_dep
                          AND det.dt < sta.dt_ret",
                      drv = "SQLite")
  cat("Done.\n")
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

  #Create spatial object
  cat("\nCreating spatial geometry... ")
  res <- res %>%
    sf::st_as_sf(agr = "identity", coords = c("x", "y"), crs = crs, remove = FALSE)
  cat("Done.\n")

  #Set the S3 class
  class(res) <- c("dets", class(res))

  #Return
  return(res)
}

#' Check if object is of class \code{dets}.
#'
#' @describeIn proc_dets Checks whether an object is of class \code{dets}.
#'
#' @param x An object to check with \code{is.dets()}.
#'
#' @export
is.dets <- function(x) {
  inherits(x, "dets")
}

#' Plots processed station history
#'
#' This function takes a processed detection history returned by
#' \code{\link{proc_dets}()} and plots the station history.
#'
#' @usage plot_sta_history(proc_det, set_par = TRUE, use_ggplot = FALSE, ...)
#'
#' @param proc_det A \code{data.frame} of class \code{dets} as returned
#' by the function \code{\link{proc_dets}()}.
#' @param set_par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? (This should be \code{FALSE} only if the
#' user wishes to manually set the graphical parameters).
#' @param use_ggplot \code{TRUE} or \code{FALSE}. Should the plot be
#' made using \code{ggplot2}? Defaults to \code{FALSE} and uses base R
#' plotting.
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
#' #Plot station history with ggplot2
#' library(ggplot2)
#' (p <- plot_sta_history(proc.det, use_ggplot = TRUE)) #ggplot
#' #Change options on ggplot
#' p +
#' theme_bw() +
#' labs(title = "Detection History")
#'
#' @export
plot_sta_history <- function(proc_det, set_par=TRUE, use_ggplot = FALSE, ...) {
  #Check class
  if(!is.dets(proc_det)){
    stop("\n  Object proc_det must be of class 'dets'.
         See ?proc_dets for conversion.")
  }

  #Factor the station id
  proc_det$sta_id <- factor(proc_det$sta_id)
  levs <- levels(proc_det$sta_id)

  #Decide whether or not to use ggplot
  if (use_ggplot){
    ##ggplot Plotting
    p <- ggplot2::ggplot(data = proc_det, aes(x = dt, y = sta_id)) +
      ggplot2::geom_point() +
      ggplot2::xlab("Date") +
      ggplot2::ylab("Station")

    return(p)

  } else {
    ##Base R Plotting

    #Store original parameters
    orig.par <- par(no.readonly = TRUE)
    #Change parameters
    if(set_par){
      par(mar = c(2.5, 6.5, 2.0, 1.0),
          cex = 0.8)
    }

    #Plot
    {plot(x = proc_det$dt, y = proc_det$sta_id, pch = 16,
          axes = FALSE, xlab=NA, ylab=NA, main="Detections by Station")
      axis.POSIXct(1, proc_det$dt)
      axis(2, labels=levs, at=1:length(levs), las=2)
      box()}
    #Revert to original parameters
    if(set_par){
      on.exit(par(orig.par), add = TRUE)
    }
  }

}

#' Plots a map of station locations
#'
#' This function takes a processed detection history returned by
#' \code{\link{proc_dets}()} and plots a map showing the relative number
#' of detections per station.
#'
#' @usage map_dets(proc_det, base_layers = NULL,
#'              base_borders = NULL, base_cols = NULL,
#'              sta_col = "black", sta_bg = "red",
#'              leg_pos = "bottomleft", set_par = TRUE, ...)
#'
#' @param proc_det A \code{data.frame} of class \code{dets} as returned
#' by the function \code{\link{proc_dets}()}.
#' Defaults to integer 4326, indicating EPSG code 4326 for WGS 84 lat/long.
#' @param base_layers A \code{list} of georeferenced vectors to be plotted
#' as the basemap under the stations. Designed to be \code{\link[sf]{sfc}}
#' \code{LINESTRING} or \code{POLYGON} objects.
#' @param base_borders A \code{list} of colors to used to outline the
#' polygons listed in \code{base_layers}. Should be the same length as
#' \code{base_layers}.
#' @param base_cols A \code{list} of colors to used to fill the polygons
#' listed in \code{base_layers}.  Should be the same length as \code{
#' base_layers}.
#' @param leg_pos Indicates where the legend should be. Can be any \code{x}
#' that the base \code{graphics} function \code{\link[graphics]{legend}}
#' will accept. Defaults to the text string \code{"bottomleft"}. Not
#' applicable if \code{use_ggplot} is \code{TRUE}.
#' @param xlim Vector of length 2 giving the limits of the plot along
#' the x-axis (longitude or easting).
#' @param ylim Vector of length 2 giving the limits of the plot along the
#' y-axis (latitude or northing).
#' @param set_par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? This should be \code{FALSE} if the
#' user wishes to manually set the graphical parameters, \emph{e.g.}, so
#' that they can still add to the plot manually after it is generated. Not
#' applicable if \code{use_ggplot} is \code{TRUE}.
#' @param return_df \code{TRUE} or \code{FALSE}. Should the function also
#' return a \code{data.frame} of aggregated station detections?
#' @param use_ggplot \code{TRUE} or \code{FALSE}. Should the plot be
#' made using \code{ggplot2}? Defaults to \code{FALSE} and uses base R
#' plotting. \emph{Note}: if both \code{return_df} and \code{use_ggplot}
#' are \code{TRUE}, then function returns a \code{list}. See details.
#' @param ... Additional arguments (not currently implemented)
#'
#' @details Details here.
#'
#' Note that if both \code{return_df} and \code{use_ggplot}
#' are \code{TRUE}, then function returns a \code{list} with elements
#' \code{df} containing the data.frame of aggregated detections and
#' \code{ggplot} containing the \code{ggplot} object for the user to
#' further manipulate.
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
#' map_dets(proc_det = proc.det,
#'     base_layers = list(acoustic$study_area, acoustic$land),
#'     base_cols = list("gray30", "wheat"), base_borders = list(NA, "seagreen"),
#'     xlim=c(-64.628, -64.612), ylim=c(17.770, 17.795),
#'     leg_pos = "topleft", sta_col = "blue", sta_bg = "orange")
#'
#' @export
map_dets <- function(proc_det, base_layers = NULL,
                     base_borders = NULL, base_cols = NULL,
                     sta_col = "black", sta_bg = "red",
                     leg_pos = "bottomleft",
                     xlim = NULL, ylim = NULL, set_par = TRUE,
                     return_df = TRUE, use_ggplot = FALSE, ...) {
  #Check class
  if(!is.dets(proc_det)){
    stop("\n  Object proc_det must be of class 'dets'.
  See ?proc_dets for conversion.")
  }

  ##Data pre-processing
  #Aggregate station detections
  stas <- proc_det %>%
    group_by(sta_id, x, y) %>%
    summarize(n_det = n()) %>%
    ungroup() %>%
    mutate(sz = cut(n_det, breaks=pretty(n_det), labels=paste0("<", pretty(n_det)[-1])))

  #Find xlim and ylim by checking extent of all layers
  #Combine all layers in a list
  if(is.null(base_layers)){
    all.layers <- list()
  } else {
    all.layers <- base_layers
  }
  all.layers[[length(base_layers)+1]] <- stas

  if(!is.null(base_layers)){
    #First check that the color lists are the right lengths
    ##Borders
    if (!is.null(base_borders)){
      #Make sure it is same length as base_layers
      if(length(base_layers) != length(base_borders)){
        stop("If you specify 'base_borders', the list must be the
               same length as 'base_layers'. Note that you can pad this
               list with NAs where necessary.")
      }
    } else {
      base_borders <- rep(NA, length(base_layers))
    }
    ##Colors
    if(!is.null(base_cols)){
      #Make sure it has the same length as base_layers
      if(length(base_layers) != length(base_cols)){
        stop("If you specify 'base_cols', the list must be the
               same length as 'base_layers'. Note that you can pad this
               list with NAs where necessary.")
      }
    } else {
      base_cols <- rep(NA, length(base_layers))
    }
  }

  if (is.null(xlim)){  #If xlim is not specified, pull from layers
    #xmin
    xmin <- min(unlist(lapply(all.layers, function(x) sf::st_bbox(x)["xmin"])))
    #xmax
    xmax <- max(unlist(lapply(all.layers, function(x) sf::st_bbox(x)["xmax"])))
    #Combine
    xlim <- c(xmin, xmax)
  }

  if (is.null(ylim)){  #If ylim is not specified, pull from layers
    #ymin
    ymin <- min(unlist(lapply(all.layers, function(x) sf::st_bbox(x)["ymin"])))
    #ymax
    ymax <- max(unlist(lapply(all.layers, function(x) sf::st_bbox(x)["ymax"])))
    #Combine
    ylim <- c(ymin, ymax)
  }

  #Decide whether to use ggplot or base
  if(use_ggplot){
    #ggplot Plotting

    #Start the plot with the station data
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = stas, aes(size = sz),
              pch = 21, color = sta_col, fill = sta_bg,
              show.legend = "point")

    #Add basemap if applicable
    if(!is.null(base_layers)){
      for(i in 1:length(base_layers)){
        #If it's a LINESTRING
        if (any(class(base_layers[[i]])=="sfc_LINESTRING")){
          p <- p +
            ggplot2::geom_sf(data = base_layers[[i]],
                    color = ifelse(is.na(base_cols[[i]]), 1, base_cols[[i]]))
        }

        #If it's a POLYGON
        if (any(class(base_layers[[i]])=="sfc_POLYGON")){
          p <-ggplot2:: p +
            geom_sf(data = base_layers[[i]],
                    color = ifelse(is.na(base_borders[[i]]), 1, base_borders[[i]]),
                    fill = ifelse(is.na(base_cols[[i]]), 1, base_cols[[i]]))
        }
      }
    }

    #Now change theme, limits, etc.
    p <- p +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.5) +
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::labs(size = "Detections") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())

    #Decide what to return
    if(return_df){
      #Must return a list, but first print plot
      ggplot2::print.ggplot(p)
      out <- list(df = stas, ggplot = p)
      return(out)
    } else {
      return(p)
    }

  } else {
    ##Base R Plotting

    #Store original parameters
    orig.par <- par(no.readonly = TRUE)
    #Change parameters
    if(set_par){
      par(mar = c(0.1, 0.1, 0.1, 0.1),
          cex = 0.8)
    }

    #Create sf point object
    dummy.point <- sf::st_point(x = c(mean(xlim), mean(ylim)), dim="XY")
    plot(dummy.point, xlim = xlim, ylim = ylim, pch=NA)


    if(!is.null(base_layers)){
      #Now that we have the colors, plot each layer
      #Loop through each layer
      for (i in 1:length(base_layers)){

        #Should be either "sfc_LINESTRING" or "sfc_POLYGON"
        if(!any(class(base_layers[[i]]) %in%
                c("sfc_LINESTRING", "sfc_POLYGON"))){

          stop(paste("Element", i, "of list 'base_layers' is neither an
                    'sfc_LINESTRING' nor 'sfc_POLYGON'."))

        } else { #As long as the classes are right

          #If it's a LINESTRING
          if (any(class(base_layers[[i]])=="sfc_LINESTRING")){
            plot(base_layers[[i]], add=TRUE,
                 col=ifelse(is.na(base_cols[[i]]), 1, base_cols[[i]]))
          }

          #If it's a POLYGON
          if (any(class(base_layers[[i]])=="sfc_POLYGON")){
            plot(base_layers[[i]], add=TRUE,
                 col=base_cols[[i]],
                 border=ifelse(is.na(base_borders[[i]]), 1, base_borders[[i]]))
          }
        }
      }
    } #End of basemap plotting

    #Plot stations
    plot(stas$geometry, add = TRUE,
         pch = 21, col = sta_col, bg = sta_bg,
         cex = as.numeric(stas$sz)/2)
    #Plot legend
    legend(x = leg_pos, legend = levels(stas$sz),
           pch = 21, col = sta_col, pt.bg = sta_bg,
           pt.cex = (1:length(levels(stas$sz)))/2,
           title = expression(bold("Detections")))

    #Revert to original parameters
    if(set_par){
      on.exit(par(orig.par), add = TRUE)
    }

    #Return the aggregated data
    return(sf::st_drop_geometry(stas))
  }




}

#' Plot a \code{dets} object
#'
#' S3 method for plotting a \code{dets} object.
#'
#' @param proc_det A \code{data.frame} of class \code{dets} as returned
#' by the function \code{\link{proc_dets}()}.
#' @param which A character string for which plot to create. Defaults
#' to \code{map}. See details below.
#' @param ... Additional parameters to pass to the specific plotting function.
#'
#' @details
#' Wrapper that calls one of the plotting functions, depending on the value
#' of parameter \code{which}. If \code{which == "history"} then the function
#' \code{\link{plot_sta_history}()} is called (the default), but if
#' \code{which == "map"} then the function \code{\link{map_dets}()} is called.
#'
#' @examples
#' #Load the example data set
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Plot station history
#' plot(proc.det)
#'
#' #Map detections
#' plot(proc.det, which = "map")
#' @export
plot.dets <- function(proc_det, which = "history", ...){
  if(!(which %in% c("map", "history"))){
    stop("\n  which must be one of 'map' (default) or 'history'")
  }
  if(which == "history"){
    plot_sta_history(proc_det, ...)
  }
  if(which == "map"){
    map_dets(proc_det, ...)
  }
}

