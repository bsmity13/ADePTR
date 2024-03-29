#Script contains descriptions of data saved in "data/"

#' Sample acoustic detection data set
#'
#' Sample files containing stations, detections, a polygon of the study area
#' boundary, and a polygon of land (barrier to movement).
#'
#' @usage \code{data(acoustic)}
#'
#' @format A list containing 4 items:
#' \describe{
#'   \item{\strong{stations}}{A \code{data.frame} relating the receivers to the stations where
#'   they were located
#'     \itemize{
#'       \item{\strong{sta_id} -- A unique identifier for each station}
#'       \item{\strong{rec_id} -- A unique identifier for each receiver}
#'       \item{\strong{dt_dep} -- The date and time that receiver was deployed
#'       at the station}
#'       \item{\strong{dt_ret} -- The date and time that receiver was
#'       retrieved from the station}
#'       \item{\strong{x} -- The x-coordinate (longitude) for the station}
#'       \item{\strong{y} -- The y-coordinate (latitude) for the station}
#'      }}
#'    \item{\strong{detections}}{A \code{data.frame} containing the acoustic
#'    detections
#'      \itemize{
#'        \item{\strong{id} -- A unique identifier for each tag deployed on
#'        a study animal}
#'        \item{\strong{dt} -- The date and time of the detection}
#'        \item{\strong{rec_id} -- The receiver at which the detection occurred}
#'      }}
#'    \item{\strong{study_area}}{A polygon of class \code{\link[sf:sfc]{sfc}}
#'    from the package \code{sf} that delimits the study area}
#'    \item{\strong{land}}{A polygon of class \code{\link[sf:sfc]{sfc}} from
#'    the package \code{sf} that designates land (a barrier to movement)}
#'   }
"acoustic"
