#Step 5
#Functions to interpolate locations and animate paths

#Conductance Raster----

#' Create a conductance raster
#'
#' Creates a conductance raster for the study area for use in calculating
#' least cost paths between two locations.
#'
#' @usage raster_conduct(sa_rast, impassable, impassable_value = 100)
#'
#' @param sa_rast A \code{\link[raster:Raster-class]{RasterLayer}} object
#' on which to calculate conductance.
#' @param impassable An \code{\link[sf]{sfc_POLYGON}} object that contains
#' land or other impassable terrain.
#' @param impassable_value A \emph{resistance} value to assign to the
#' \code{impassable} cells of the raster. This should be a large number
#' to ensure the animal will not cross the impassable sections. The function
#' will convert it from a resistance to a conductance.
#'
#' @return Returns a \code{\link[raster:Raster-class]{RasterLayer}} object
#' of conductance values.
#'
#' @details This function takes a \code{\link[raster:Raster-class]{RasterLayer}}
#' as a template and an \code{\link[sf]{sfc_POLYGON}} and returns another
#' \code{RasterLayer} of conductance. This conductance raster is intended to
#' be passed downstream for calculation of least cost paths between animal
#' locations.
#'
#' @examples
#' #Load sample data
#' data(acoustic)
#'
#' #Initialize raster
#' r.cond <- init_raster(study_area = acoustic$study_area, res = 0.001)
#'
#' #Create conductance raster
#' conduct <- raster_conduct(sa_rast = r.cond, impassable = acoustic$land)
#'
#' #Plot
#' plot(conductance)
#' plot(acoustic$study_area, add = TRUE)
#' plot(acoustic$land, add = TRUE, col = NA)
#'
#' @export
raster_conduct <- function(sa_rast, impassable, impassable_value = 100){

  #Convert impassable from sf to sp
  imp_sp <- as(impassable, "Spatial")

  #Rasterize resistance
  resist <- raster::rasterize(x = imp_sp, y = sa_rast,
                              fun = function(x, na.rm){
                                return(impassable_value)},
                              background = 1)
  #Convert to conductance
  conduct <- 1/resist

  return(conduct)
}

#Linear paths----

#' Connect successive locations with lines
#'
#' Creates straight path line objects between each location for
#' any number of individuals.
#'
#' @export
str_path <- function(det_locs){

  #Make sure "id" is in the data.frame of det_locs
  if(!("id" %in% names(det_locs))){
    stop("There must be a column named 'id' that identifies individuals.
           See ?proc_dets for proper formatting.")
  }

  paths <- det_locs %>%
    arrange(id, dt) %>%
    group_by(id) %>%
    rename(x1 = x, y1 = y) %>%
    mutate(x2 = lead(x1, n = 1L),
           y2 = lead(y1, n = 1L)) %>%
    filter(!is.na(x2)) %>%
    rowwise() %>%
    mutate(geometry = sf::st_geometry(
      sf::st_linestring(c(
        sf::st_point(c(x1, y1)),
        sf::st_point(c(x2, y2))
        ))))

  return(paths)
}

#Least Cost Paths----

#' Create least cost paths between each location
#'
#' Creates least cost paths as line objects between each location for
#' a single individual.
#'
#' @param det_locs_ind A \code{data.frame} of georeferenced acoustic
#' detections, subset for a single individual. Could be, \emph{e.g.},
#' the \code{data.frame} returned by \code{\link{proc_dets}()} or the
#' \code{data.frame} returned by \code{\link{coa_locs}()}, subset to a
#' single individual.
#' @param tl A \code{\link[gdistance:transition]{TransitionLayer}}
#' object used to calculate least cost path.
#'
#' @return Returns an boject of class \code{lcp_ind}, which is also
#' a \code{tbl_df}, \code{tbl}, and \code{data.frame}.
#'
#' @details This function is intended to be called internally by the
#' function \code{\link{lc_path}()}.
#'
#' @export
lc_path_ind <- function(det_locs_ind, tl){

  cat(paste0("\n  Calculating path for ID: ", unique(det_locs_ind$id)))

  lcp_ind <- det_locs_ind %>%
    arrange(dt) %>%
    mutate(row = row_number()) %>%
    filter(!is.na(x2)) %>%
    group_by(row) %>%
    mutate(geom = sf::st_as_sfc(
      gdistance::shortestPath(tl,
                              c(x1, y1),
                              c(x2, y2),
                              output = "SpatialLines"))) %>%
    ungroup() %>%
    mutate(n_pts = sapply(geom, length)/2)

  class(lcp_ind) <- c("lcp_ind", class(lcp_ind))

  return(lcp_ind)
}

#' Create least cost paths between each location
#'
#' Creates least cost paths as line objects between each location for
#' multiple individuals.
#'
#' @param det_locs A \code{data.frame} of georeferenced acoustic
#' detections, subset for a single individual. Could be, \emph{e.g.},
#' the \code{data.frame} returned by \code{\link{proc_dets}()} or the
#' \code{data.frame} returned by \code{\link{coa_locs}()}, subset to a
#' single individual.
#' @param conductance A \code{\link[raster:Raster-class]{RasterLayer}}
#' object of conductance, either returned by \code{\link{raster_conduct}()}
#' or a similar object created by the user.
#'
#' @return Returns an object of class \code{lcp_list}, which is also
#' a \code{list}. List elements are of class \code{lcp_ind}. See
#' \code{\link{lc_path_ind}} for details.
#'
#' @details Estimating the least cost paths between a large number of
#' locations can be time consuming. To help with this, the function
#' reports to the user its progress by printing to the console which
#' ID it is currently working on.
#'
#' @examples
#'
#' #Load sample data
#' data(acoustic)
#'
#' #Process detections and stations
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Calculate  COAs -- Note, using 6 hour COAs to reduce processing time
#' coas.6h <- coa_locs(proc.det, Delta_t = "6 hours")
#'
#' #Initialize raster
#' r.cond <- init_raster(study_area = acoustic$study_area, res = 0.001)
#'
#' #Create conductance raster
#' conduct <- raster_conduct(sa_rast = r.cond, impassable = acoustic$land)
#'
#' #Create least-cost paths
#' lcps <- lc_path(det_locs = coas.6h, conductance = conduct)
#'
#' @export
lc_path <- function(det_locs, conductance){

  #Make sure "id" is in the data.frame of det_locs
  if(!("id" %in% names(det_locs))){
    stop("There must be a column named 'id' that identifies individuals.
           See ?proc_dets for proper formatting.")
  }

  #Create transition layer
  tl <- gdistance::transition(x = conductance,
                              transitionFunction = mean,
                              directions = 8)
  #Apply geoCorrection
  tl <- geoCorrection(tl, type = "c", multpl = FALSE)

  #Count how many unique ids there are
  n_ids <- length(unique(det_locs$id))

  #Report to user
  cat(paste0("\nConstructing least-cost paths for ", n_ids, " IDs..."))

  #Process
  paths <- det_locs %>%
    arrange(id, dt) %>%
    group_by(id) %>%
    rename(x1 = x, y1 = y) %>%
    mutate(x2 = lead(x1, n = 1L),
           y2 = lead(y1, n = 1L)) %>%
    filter(!is.na(x2)) %>%
    split(.$id) %>%
    purrr::map(lc_path_ind, tl = tl)

  class(paths) <- c("lcp_list", "list")

  return(paths)
}

#' Remove length-0 lines from an \code{lcp_ind} object
#'
#' Takes an object of class \code{lcp_ind} and filters out any
#' lines that do not have at least 2 vertices.
#'
#' @usage filter_line0(lcp_ind)
#'
#' @param lcp_ind An object of class \code{lcp_ind}.
#'
#' @return Returns an object of class \code{lcp_ind}.
#'
#' @details This function is intended for internal use on \code{lcp_list}
#' objects.
#'
#' @export
filter_line0 <- function(lcp_ind){
  lcp_ind %>%
    filter(n_pts > 1) %>%
    return()
}

#' Plot an object of class \code{lcp_list}
#'
#' Plots the least cost paths between locations for all individuals
#' in an \code{lcp_list} object.
#'
#' @usage plot_lc_path(...)
#'
#' @param lcps An object of class \code{lcp_list}
#' @param path_palette A function name that will generate a palette for each
#' of the individuals tracked. Can be any function that will accept \code{n}
#' as an argument.
#' @param add \code{TRUE} or \code{FALSE}. Should the plot be added to an
#' existing plot? Defaults to \code{FALSE}. Not applicable if
#' \code{use_ggplot = TRUE}. User could, \emph{e.g.}, use this option to
#' add this plot on top of \code{\link{map_dets}()}.
#' @param use_ggplot \code{TRUE} or \code{FALSE}. Should the plot be
#' made using \code{ggplot2}? Defaults to \code{FALSE} and uses base R
#' plotting.
#'
#' @details Blah blah
#'
#' @examples
#'
#' #Load sample data
#' data(acoustic)
#'
#' #Process detections and stations
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Calculate  COAs -- Note, using 6 hour COAs to reduce processing time
#' coas.6h <- coa_locs(proc.det, Delta_t = "6 hours")
#'
#' #Initialize raster
#' r.cond <- init_raster(study_area = acoustic$study_area, res = 0.001)
#'
#' #Create conductance raster
#' conduct <- raster_conduct(sa_rast = r.cond, impassable = acoustic$land)
#'
#' #Create least-cost paths
#' lcps <- lc_path(det_locs = coas.6h, conductance = conduct)
#'
#' #Plot
#' plot_lc_path(lcps)
#'
#' #Plot on top of detections
#'   #Note set_par = FALSE in both functions
#'   #Consider manually setting graphical parameters for this plot
#'
#' map_dets(proc_det = proc.det, sta_crs = 4326,
#'   base_layers = list(acoustic$study_area, acoustic$land),
#'   base_cols = list("gray30", "wheat"), base_borders = list(NA, "seagreen"),
#'   xlim=c(-64.636, -64.604), ylim=c(17.770, 17.795),
#'   leg_pos = "topleft", sta_col = "blue", sta_bg = "orange",
#'   set_par = FALSE)
#'
#' plot_lc_path(lcps, set_par = FALSE, add = TRUE)
#'
#' #Plot with ggplot
#' plot_lc_path(lcps, use_ggplot = TRUE) +
#'   theme_bw() +
#'   geom_sf(data = acoustic$land)


#' @export
plot_lc_path <- function(lcps, path_palette = viridis::viridis, add = FALSE, set_par = !add,
                         use_ggplot = FALSE){
  #Check class
  if(!("lcp_list" %in% class(lcps))){
    stop("Object must be of class 'lcp_list'.\n  See ?lc_path for details.")
  }

  #Remove 0-length lines
  lcps <- lcps %>%
    purrr::map(filter_line0)

  #Process colors
  path_cols <- do.call(path_palette, args = list(n = length(lcps)))

  #Decide whether to use ggplot or base plotting
  if(use_ggplot){
    #ggplot Plotting
    #Convert to data.frame
    suppressWarnings(lcps_df <- bind_rows(lcps))
    #Reformat geom as sf object
    lcps_df$geometry <- sf::st_as_sfc(lcps_df$geom)
    ggp <- ggplot2::ggplot(data = lcps_df) +
      ggplot2::geom_sf(aes(color = id)) +
      scale_color_manual(values = path_cols)
    return(ggp)
  } else {
    #Base R plotting
    #Store original parameters
    orig.par <- par(no.readonly = TRUE)
    #Change parameters
    if(set_par){
      par(mar = c(0.1, 0.1, 0.1, 0.1),
          cex = 0.8)
    }
    #Loop through the lcps
    for (i in 1:length(lcps)){
      if (i == 1){
        plot(lcps[[i]]$geom, col = path_cols[[i]], add = add)
      } else {
        plot(lcps[[i]]$geom, col = path_cols[[i]], add = TRUE)
      }
    }
    #Revert to original parameters
    if(set_par){
      on.exit(par(orig.par), add = TRUE)
    }
  }
}

#' Plot an object of class \code{lcp_list}
#'
#' Plots the least cost paths between locations for all individuals
#' in an \code{lcp_list} object.
#'
#' @details See \code{\link{plot_lc_path}()} for details.
#'
#' @examples
#'
#' #Load sample data
#' data(acoustic)
#'
#' #Process detections and stations
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Calculate  COAs -- Note, using 6 hour COAs to reduce processing time
#' coas.6h <- coa_locs(proc.det, Delta_t = "6 hours")
#'
#' #Initialize raster
#' r.cond <- init_raster(study_area = acoustic$study_area, res = 0.001)
#'
#' #Create conductance raster
#' conduct <- raster_conduct(sa_rast = r.cond, impassable = acoustic$land)
#'
#' #Create least-cost paths
#' lcps <- lc_path(det_locs = coas.6h, conductance = conduct)
#'
#' #Generic plot
#' plot(lcps)
#'
#' #Generic method with land added
#' plot(lcps, set_par = FALSE)
#' plot(acoustic$land, add = TRUE, col = "wheat")
#'
#' @export
plot.lcp_list <- function(lcps, ...){
  plot_lc_path(lcps, ...)
}

#Location Interpolation----






