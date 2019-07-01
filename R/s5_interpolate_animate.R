#Step 5
#Functions to interpolate locations and animate paths

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
#' r.sa <- init_raster(study_area = acoustic$study_area)
#'
#' #Create conductance raster
#' conduct <- raster_conduct(sa_rast = r.sa, impassable = acoustic$land)
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

#Calculate least cost path----


