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
    mutate(n_pts = 2) %>%
    ungroup()

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
  tl <- gdistance::geoCorrection(tl, type = "c", multpl = FALSE)

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

#' Check if object is of class \code{*_paths}.
#'
#' Convenience function that checks is object is of class \code{*_paths}.
#' Currently must be either \code{str_paths} or \code{lc_paths}.
#' @export
is.paths <- function(x) {
  return(inherits(x, "str_paths") | inherits(x, "lc_paths"))
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
  if(!is.paths(paths)){
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
      ggplot2::geom_sf(aes(geometry = geometry, color = id), show.legend = "point") +
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

#' Generate regular points along a path
#'
#' Takes a *_paths object and generates regularly-spaced
#' points along the path.
#'
#' @usage regular_points(paths, Delta_t)
#'
#' @param paths An object of class \code{*_paths}: either \code{str_paths}
#' or \code{lc_paths}
#' @param Delta_t The timestep to use for interpolation. The number of
#' interpolated points will be approximately (total time)/Delta_t. See
#' details below.
#'
#' @return Returns an object with S3 Class \code{reg_points}. Object is also
#' of class \code{sf}, \code{tbl_df}, \code{tbl}, and \code{data.frame}.
#'
#' @details This function returns points that are regularly placed in
#' \emph{space}, even though the user provides a \emph{time} argument
#' (\code{Delta_t}). When the input locations are regularly spaced in time,
#' the output locations will be regularly spaced in both
#' \emph{space and time}. This will already be the case if the user has
#' decided to use COAs to represent the locations.
#'
#' The number of points returned is reported in the output \code{data.frame}
#' in the column \code{$n_interp}. We determine this number by counting the
#' length of the sequence from the start time to the end time of the
#' track for each individual, with the \code{by} argument being
#' \code{Delta_t}. \emph{I.e.}:
#'
#' \preformatted{
#'   n_interp = length(
#'                seq(
#'                  from = min(dt),
#'                  to = max(dt),
#'                  by = Delta_t))}
#'
#' Therefore, \code{Delta_t} can be any value that can be accepted by
#' the \code{by} argument of \code{\link{seq.POSIXt}()}.
#'
#' @export
regular_points <- function(paths, Delta_t){
  #Check class
  if(!is.paths(paths)){
    stop("Object must be of class 'str_paths' or 'lc_paths'.
           See ?str_paths or ?lc_paths for details.")
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

  #Interpolate along the path
  path_summ <- paths %>%
    filter(n_pts > 1) %>%
    group_by(id) %>%
    summarize(start_dt = min(dt, na.rm = TRUE),
              end_dt = max(dt, na.rm = TRUE),
              lines = sf::st_combine(geometry)) %>%
    rowwise() %>%
    mutate(n_interp = length(
      seq(from = start_dt,
          to = end_dt,
          by = Delta_t)))

  res <- path_summ %>%
    mutate(pts = sf::st_combine(
      sf::st_sample(lines,
                    size = n_interp,
                    type = "regular"))) %>%
    dplyr::select(id, start_dt, end_dt, n_interp, geometry = pts) %>%
    ungroup()

  #Now convert from wide (MULTIPOINT) to long (POINT) data
  res <- sf::st_sf(res)
  res <- sf::st_cast(x = res, to = "POINT", warn = FALSE)

  #Add vector of datetimes to res
  res$dt <- do.call("c", lapply(path_summ$id, function(x){
    sdt <- path_summ %>%
      dplyr::filter(id == x) %>%
      dplyr::pull(start_dt)
    edt <- path_summ %>%
      dplyr::filter(id == x) %>%
      dplyr::pull(end_dt)
    npt <- path_summ %>%
      dplyr::filter(id == x) %>%
      dplyr::pull(n_interp)
    dt_seq <- seq(from = sdt, to = edt, length.out = npt)
    return(dt_seq)
  }))

  #Select columns of interest
  res <- res %>%
    dplyr::select(id, dt, geometry)

  #Set CRS to the projected CRS
  sf::st_crs(res) <- proj_crs

  #Convert to orginal CRS
  if(sf::st_crs(res) != orig_crs){
    res <- sf::st_transform(res, orig_crs)
  }

  #Assign S3 Class
  class(res) <- c("reg_points", class(res))

  return(res)
}

#' Simple animation of regular points
#'
#' Provides a simple animation of a \code{reg_points} object.
#'
#' @usage animate_points(reg_points, ani.width = 800, ani.height = 600,
#' ...)
#'
#' @param reg_points An object of class \code{reg_points} to create
#' the animation.
#' @param ani.width Width of the animation in pixels. Passed to
#' \code{\link[animation:ani.options]{animation::ani.options}()}.
#' @param ani.height Height of the animation in pixels. Passed to
#' \code{\link[animation:ani.options]{animation::ani.options}()}.
#'
#' @return Returns an object of class \code{gganim} from package
#' \link[gganimate:gganimate-package]{gganimate}.
#'
#' @details This function quickly creates a simple animation of a
#' \code{reg_points} object returned by \code{\link{regular_points}()}.
#' It uses the package \link[gganimate:gganimate-package]{gganimate}
#' to make an animation. The function returns an object of class
#' \code{gganim}, so the user can manipulate the object prior to printing
#' just as they would manipulate any ggplot with \code{`+`}. Printing the
#' object automatically renders the animation using
#' \code{\link[gganimate]{animate}}.
#'
#' We encourage the user to explore
#' \link[gganimate:gganimate-package]{gganimate} for producing their own
#' custom animations from the interpolated regular points returned by
#' \code{\link{regular_points}()}. The source code from this function
#' could serve as a starting point.
#'
#' \preformatted{
#'     ani <- ggplot(reg_points) +
#'       geom_sf(aes(color = id), key_glyph = "point") +
#'       theme_bw() +
#'       transition_time(time = dt) +
#'       labs(title = 'Time: {frame_time}', color = "ID")
#' }
#'
#' @export
animate_points <- function(reg_points,
                           ani.width = 800, ani.height = 600, ...){
  #Check class
  if(!inherits(reg_points, "reg_points")){
    stop("Object must be of class 'reg_points'.
           See ?regular_points for details.")
  }

  #Set animation options
  animation::ani.options(ani.width = ani.width, ani.height = ani.height)

  #Define animation
  ani <- ggplot(reg_points) +
    geom_sf(aes(color = id), key_glyph = "point") +
    theme_bw() +
    transition_time(time = dt) +
    labs(title = 'Time: {frame_time}', color = "ID")

  #Return
  return(ani)
}



