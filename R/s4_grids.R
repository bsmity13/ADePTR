#Step 4
#Summarize detections on a grid

#' Creates an empty raster for detection summaries
#'
#' This function can take one of user-defined x- and y-limits, station
#' locations, or a study area polygon and create an empty raster to use
#' in the detection summary figures.
#'
#' @usage init_raster(xlim = NULL, ylim = NULL, sta = NULL,
#'                    study_area = NULL, res = 0.01,
#'                    crs = sp::CRS("+init=epsg:4326"))
#'
#' @param xlim A numeric vector of length 2 giving user-defined x-limits. If
#' \code{xlim} is provided, \code{ylim} must also be provided. These fields
#' override both \code{sta} and \code{study_area}.
#' @param ylim A numeric vector of length 2 giving user-defined y-limits. If
#' \code{xlim} is provided, \code{ylim} must also be provided. These fields
#' override both \code{sta} and \code{study_area}.
#' @param sta A data.frame with x- and y-coordinates of all the stations
#' in the study area. The x- and y-coordinates should be in columns named
#' \code{"x"} and \code{"y"}, respectively. Ignored if \code{xlim} or
#' \code{ylim} are provided. This field overrides \code{study_area}.
#' @param study_area A polygon of class \code{\link[sf:sfc]{sfc}}
#' from the package \code{sf} that delimits the study area.
#' @param res A numeric vector  of length 1 or 2 giving the size of a grid
#' cell in units of the coordinate reference system. See
#' \code{\link[raster]{res}} for details.
#' @param crs A character or object of class CRS. Defines the coordinate
#' reference system to use for the resulting raster. Should match the CRS of
#' the input if applicable. Used by \code{\link[raster]{raster}} as
#' the \code{crs} argument. Defaults to \code{(CRS("+init=epsg:4326"))},
#' indicating EPSG code 4326 for WGS 84 lat/long.
#'
#' @details Behavior depends upon the input. Once the function reaches a valid
#' input, it returns the raster without bothering to check if additional inputs
#' were provided.
#'
#' The order it checks in is:
#' \enumerate{
#'   \item{\code{xlim} and \code{ylim}} (must provide both)
#'   \item{\code{sta}} (x- and y-coordinates should be named \code{"x"} and
#'   \code{"y"})
#'   \item{\code{study_area}} (polygon of class \code{\link[sf:sfc]{sfc}})
#' }
#'
#' @return Returns an object of class
#' \code{\link[raster:Raster-class]{RasterLayer}}
#'
#' @examples
#' #Load sample data
#' data(acoustic)
#'
#' #Using xlim and ylim
#' xlim <- c(-64.628, -64.612)
#' ylim <- c(17.770, 17.795)
#' r.lim <- init_raster(xlim = xlim, ylim = ylim)
#'
#' #Use UTM
#' xlim <- c(321160, 335010)
#' ylim <- c(1965700, 1973810)
#' r.lim.utm <- init_raster(xlim = xlim, ylim = ylim,
#'                          res = 500, #units here are meters
#'                          crs = sp::CRS("+init=epsg:32620"))
#'
#' #Using sta
#' r.sta <- init_raster(sta = acoustic$stations)
#'
#' #Using study_area
#' r.sa <- init_raster(study_area = acoustic$study_area)
#'
#' #Change the resolution (in CRS units)
#' r.sa2 <- init_raster(study_area = acoustic$study_area,
#'                      res = 0.0001)
#'
#'
#'@export
init_raster <- function(xlim = NULL, ylim = NULL, sta = NULL,
                        study_area = NULL, res = 0.005,
                        crs = sp::CRS("+init=epsg:4326")) {

  #Make sure at least one of the inputs was provided.
  if(all(is.null(c(sta, study_area, xlim, ylim)))){
    stop("You must provide at least one of the inputs to create the raster.")
  }

  #Check each input in order to create raster

  #1. xlim/ylim
  if(!(is.null(xlim) & is.null(ylim))){

    #Make sure both are provided
    if(is.null(xlim) | is.null(ylim)){
      stop("You must provide both 'xlim' and 'ylim' if you choose that option.")
    }

    #Create extent using xlim and ylim
    ext <- raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])

  } else{
    #2. sta
    if(!is.null(sta)){
      #Get limits from sta
      xlim <- c(min(sta$x, na.rm = TRUE), max(sta$x, na.rm = TRUE))
      ylim <- c(min(sta$y, na.rm = TRUE), max(sta$y, na.rm = TRUE))

      #Create extent using xlim and ylim
      ext <- raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])

    } else {
      #3. study_area
      if(!is.null(study_area)){
        #Get the bounding box
        bb <- sf::st_bbox(study_area)

        #Create extent
        ext <- raster::extent(bb[["xmin"]], bb[["xmax"]],
                              bb[["ymin"]], bb[["ymax"]])
      }
    }
  }

  #Rasterize the extent
  rast <- raster::raster(x = ext, resolution = res, crs = crs)

  #Return
  return(rast)
}

#' Rasterize acoustic detections
#'
#' This function takes either detections at stations or centers of
#' activity (COAs) and creates a raster map showing number of detections
#' in each cell.
#'
#' @usage rasterize_dets(det_locs, r = init_raster(sta = det_locs))
#'
#' @param det_locs A \code{data.frame} of georeferenced acoustic detections.
#' Could be, \emph{e.g.}, the \code{data.frame} returned by
#' \code{\link{proc_dets}()} or the \code{data.frame} returned by
#' \code{\link{coa_locs}()}
#' @param r A blank raster to use as the template for rasterizing the
#' locations. Defaults to a call to \code{\link{init_rast}()} with
#' \code{det_locs} passed as the argument \code{sta}
#'
#' @details Calls the function \code{\link{raster::rasterize}()} to generate
#' the raster. Uses the function \code{'count'}
#'
#' @return Returns an object of class
#' \code{\link[raster:Raster-class]{RasterLayer}} with values equal
#' to the number of detections in each raster cell.
#'
#' @examples
#' #Load example data
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Default rasterize
#' r1 <- rasterize_dets(det_locs = proc.det)
#'
#' #Plot
#' plot(r1)
#' plot(acoustic$land, add = TRUE, col="wheat")
#'
#' #Rasterize with a smaller grid
#' r2 <- rasterize_dets(det_locs = proc.det,
#'                      r = init_raster(sta = proc.det, res = 0.001))
#'
#' #Plot
#' plot(r2)
#' plot(acoustic$land, add = TRUE, col="wheat")
#'
#' #Rasterize COAs
#' coas <- coa_locs(proc.det)
#' r3 <- rasterize_dets(det_locs = coas,
#'                      r = init_raster(sta = proc.det, res = 0.001))
#'
#' #Plot
#' plot(r3)
#' plot(acoustic$land, add = TRUE, col="wheat")
#'
#' @export
rasterize_dets <- function(det_locs, r = init_raster(sta = det_locs)){
  #Internal note: this function needs to be updated so that it checks
    #any 'r' provided by the user

  #Drop geometry
  det_locs <- sf::st_drop_geometry(det_locs)

  #Grab x and y coords of locations and coerce to matrix
  locs <- as.matrix(det_locs[,c("x", "y")])

  out.r <- raster::rasterize(x = locs,
                     y = r,
                     fun = "count",
                     background = NA)
  return(out.r)
}

#' Create a presence raster
#'
#' Converts a raster of counts to a raster of presence-only or
#' presence-absence
#'
#' @usage raster_pres(r)
#'
#' @param r A raster
#'
#' @details This function sets any values >= 1 to 1. If the background
#' values of the raster are NA, they remain NA in the output, allowing
#' the raster to be interpreted as presence-only. If the background
#' values of the raster are 0, then the result will be also contain 0s,
#' allowing the raster to be interpreted as presence-abscence.
#'
#' @examples
#'
#' #Load example data
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Create COAs
#' coas.60 <- coa_locs(proc.det)
#'
#' #Initialize a base raster
#' r.sa <- init_raster(study_area = acoustic$study_area, res = 0.005)
#'
#' #Rasterize COAs
#' coa.r <- rasterize_dets(det_locs = coas.60, r = r.sa)
#'
#' #Presence-only layer for all animals in the dataset
#' coa.pres <- raster_pres(coa.r)
#' #Plot
#' plot(coa.pres)
#'
#' ##Use tidy workflow to create presence for multiple individuals
#'
#' library(dplyr)
#'
#' #Animal 1
#' ani1 <- coas.60 %>%
#' filter(id == unique(id)[1]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 2
#' ani2 <- coas.60 %>%
#' filter(id == unique(id)[2]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 3
#' ani3 <- coas.60 %>%
#' filter(id==unique(id)[3]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 4
#' ani4 <- coas.60 %>%
#' filter(id==unique(id)[4]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' @export
raster_pres <- function(r){
  #Get values
  vals <- raster::getValues(r)
  #Convert to presence only
  vals[vals > 0] <- 1
  #Set values
  r.out <- raster::setValues(r, vals)
}


#' Create an overlap raster
#'
#' Takes multiple rasters of presence-only or presence-absence and
#' calculates the overlap.
#'
#' @usage raster_overlap(..., background = NA)
#'
#' @param ... Any number of \code{\link[raster:Raster-class]{RasterLayer}}
#' objects to sum.
#' @param background A value to replace any 0 values with. Either \code{0}
#' or \code{NA} (the default) are the most logical choices.
#'
#' @details This function sums the values of all the input rasters
#' and then replaces any 0 values with the background value. Note that
#' this function will accept rasters with values > 1. If that is the
#' case, the user should take caution not to interpret the values as
#' the number of overlapping individuals, but rather, as the sum of all
#' input rasters.
#'
#' @examples
#'
#' #Load example data
#' data(acoustic)
#'
#' #Process detections
#' proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)
#'
#' #Create COAs
#' coas.60 <- coa_locs(proc.det)
#'
#' #Initialize a base raster
#' r.sa <- init_raster(study_area = acoustic$study_area, res = 0.005)
#'
#' #Use tidy workflow to create presence for multiple individuals
#'
#' library(dplyr)
#'
#' #Animal 1
#' ani1 <- coas.60 %>%
#' filter(id == unique(id)[1]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 2
#' ani2 <- coas.60 %>%
#' filter(id == unique(id)[2]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 3
#' ani3 <- coas.60 %>%
#' filter(id==unique(id)[3]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Animal 4
#' ani4 <- coas.60 %>%
#' filter(id==unique(id)[4]) %>%
#' rasterize_dets(r.sa) %>%
#' raster_pres()
#'
#' #Calculate overlap
#' #Presence-only
#' ani_over <- raster_overlap(ani1, ani2, ani3, ani4)
#' plot(ani_over)
#' plot(acoustic$study_area, add=TRUE)
#' plot(acoustic$land, add=TRUE, col="wheat")
#'
#' #Presence-abscence
#' ani_over0 <- raster_overlap(ani1, ani2, ani3, ani4, background = 0)
#' plot(ani_over0)
#' plot(acoustic$study_area, add=TRUE)
#' plot(acoustic$land, add=TRUE, col="wheat")
#'
#' @export
raster_overlap <- function(..., background = NA){

  #Create list of rasters
  rs <- list(...)

  #Check that all inputs are rasters
  r_check <- unlist(lapply(rs, function(x){inherits(x, "RasterLayer")}))
  if(any(!r_check)){
    stop("All inputs must be of class 'RasterLayer'.")
  }

  #Stack rasters
  r_stack <- raster::stack(...)

  #Sum rasters
  r_sum <- sum(r_stack, na.rm = TRUE)

  #Change out background
  r_sum_vals <- getValues(r_sum)
  r_sum_vals[r_sum_vals==0] <- background
  r_sum_out <- setValues(r_sum, r_sum_vals)

  return(r_sum_out)
}
