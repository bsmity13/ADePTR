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
str_paths <- function(det_locs){

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
        )))) %>%
    mutate(n_pts = 2)

  class(paths) <- c("str_paths", class(paths))

  return(paths)
}

#Least Cost Paths----

#' Create least cost paths between each location
#'
#' Creates least cost paths as line objects between each location for
#' any number of individuals.
#'
#' @param det_locs A \code{data.frame} of georeferenced acoustic
#' detections. Could be, \emph{e.g.}, the \code{data.frame} returned by
#' \code{\link{proc_dets}()} or the \code{data.frame} returned by
#' \code{\link{coa_locs}()}.
#' @param conductance A \code{\link[raster:Raster-class]{RasterLayer}}
#' object of conductance, either returned by \code{\link{raster_conduct}()}
#' or a similar object created by the user.
#'
#' @return Returns an object of class \code{lcp_list}, which is also
#' a \code{list}. List elements are of class \code{lcp_ind}. See
#' \code{\link{lc_path_ind}} for details.
#'
#' @details Estimating the least cost paths between a large number of
#' locations can be time consuming.
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
#' lcps <- lc_paths(det_locs = coas.6h, conductance = conduct)
#'
#' @export
lc_paths <- function(det_locs, conductance){

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

  paths <- det_locs %>%
    arrange(id, dt) %>%
    group_by(id) %>%
    rename(x1 = x, y1 = y) %>%
    mutate(x2 = lead(x1, n = 1L),
           y2 = lead(y1, n = 1L)) %>%
    filter(!is.na(x2)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(geometry = sf::st_as_sfc(
      gdistance::shortestPath(tl,
                              c(x1, y1),
                              c(x2, y2),
                              output = "SpatialLines"))) %>%
    ungroup() %>%
    mutate(n_pts = sapply(geometry, length)/2)

  class(paths) <- c("lc_paths", class(paths))

  return(paths)
}

#' Plot an object of class \code{*_paths}
#'
#' Plots the paths between locations for all individuals
#' in an \code{str_paths} or an \code{lc_paths} object.
#'
#' @usage plot_paths(paths, path_palette = viridis::viridis, add = FALSE,
#' set_par = !add, use_ggplot = FALSE)
#'
#' @param paths An object of class \code{*_paths}: either \code{str_paths}
#' or \code{lc_paths}
#' @param path_palette A function name that will generate a palette for each
#' of the individuals tracked. Can be any function that will accept \code{n}
#' as an argument.
#' @param add \code{TRUE} or \code{FALSE}. Should the plot be added to an
#' existing plot? Defaults to \code{FALSE}. Not applicable if
#' \code{use_ggplot = TRUE}. User could, \emph{e.g.}, use this option to
#' add this plot on top of \code{\link{map_dets}()}.
#' @param set_par \code{TRUE} or \code{FALSE}. Should the function change
#' the graphical parameters or not? This should be \code{FALSE} if the
#' user wishes to manually set the graphical parameters, \emph{e.g.}, so
#' that they can still add to the plot manually after it is generated.
#' Defaults to \code{!add}.
#' @param use_ggplot \code{TRUE} or \code{FALSE}. Should the plot be
#' made using \code{ggplot2}? Defaults to \code{FALSE} and uses base R
#' plotting.
#' @param stps An object of class \code{str_paths} to be passed from the
#' generic \code{plot()} to \code{plot_paths()}
#' @param lcps An object of class \code{lc_paths} to be passed from the
#' generic \code{plot()} to \code{plot_paths()}
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
#' lcps <- lc_paths(det_locs = coas.6h, conductance = conduct)
#'
#' #Plot
#' plot_paths(lcps)
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
#' plot_paths(lcps, set_par = FALSE, add = TRUE)
#'
#' #Plot with ggplot
#' plot_paths(lcps, use_ggplot = TRUE) +
#'   theme_bw() +
#'   geom_sf(data = acoustic$land)


#' @export
plot_paths <- function(paths, path_palette = viridis::viridis, add = FALSE, set_par = !add,
                         use_ggplot = FALSE){
  #Check class
  if(!(any(grepl("paths", class(paths))))){
    stop("Object must be of class 'str_paths' or 'lc_paths'.
           See ?str_paths or ?lc_paths for details.")
  }

  #Remove 0-length lines
  paths <- paths %>%
    filter(n_pts > 1)

  #Process colors
  path_cols <- do.call(path_palette, args = list(n = length(unique(paths$id))))

  #Decide whether to use ggplot or base plotting
  if(use_ggplot){
    #ggplot Plotting
    ggp <- ggplot2::ggplot(data = paths) +
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
    #Setup palette for individuals
    orig.pal <- palette() #Get previous palette
    n.ind <- length(unique(paths$id))
    palette(do.call(path_palette, args = list(n = n.ind))) #Change palette
    #Plot the paths
    plot(paths$geometry,
         col = as.numeric(factor(paths$id)),
         add = add)

    #Revert to original parameters
    palette(orig.pal)
    if(set_par){
      on.exit(par(orig.par), add = TRUE)
    }
  }
}

#' Plot an object of class \code{str_paths}
#' @describeIn plot_paths Plot \code{str_paths} objects
#' @export
plot.str_paths <- function(stps, ...){
  plot_paths(stps, ...)
}

#' Plot an object of class \code{lc_paths}
#' @describeIn plot_paths Plot \code{lc_paths} objects
#' @export
plot.lc_paths <- function(lcps, ...){
  plot_paths(lcps, ...)
}

#Location Interpolation----






