#' ADePTR: Acoustic Detection Processing and Visualization Tool in R
#'
#' ADePTR is package developed to streamline processing and visualization
#' of animal movement data from a passive acoustic array.
#'
#' @docType package
#' @name ADePTR
#'
#' @details
#'
#' This package was designed with the following workflow in mind:
#'
#' \enumerate{
#'   \item Process input files
#'   \item Filter erroneous detections
#'   \item Calculate centers of activity
#'   \item Summarize detections on a grid of the study area
#' }
#'
#' For each step of the workflow, there are both possible inputs and
#' possible outputs. The package is designed to work well with minimal
#' user input (\emph{i.e.}, most arguments have defaults), but it is also
#' designed to allow a high level of customization. Advanced users should
#' be able to edit the outputs of any step before passsing it to the next
#' step, or alternatively, to enter the workflow with their own data at
#' any point.
#'
#' @section Step 1. Process Inputs:
#'
#' The functions available for this step are:
#' \itemize{
#'   \item \code{\link{proc_dets}()} -- uses station/receiver information to
#'   georeference the detections
#'   \item \code{\link{plot_sta_history}()} -- makes a summary plot of
#'   detections over time by station
#' }
#'
#' @section Step 2. Filter:
#'
#' The functions available for this step are:
#' \itemize{
#'   \item \code{\link{NOT_YET}()} -- DOES SOMETHING
#' }
#'
#'  @section Step 3. COAs:
#'
#' The functions available for this step are:
#' \itemize{
#'   \item \code{\link{coa_locs}()} -- calculates centers of activity (COAs)
#'   using either arithmetic or harmonic means
#' }
#'
#' @section Step 4. Summarize on Grid:
#'
#' The functions available for this step are:
#' \itemize{
#'   \item \code{\link{NOT_YET}()} -- DOES SOMETHING
#' }
#'
#'
#' @encoding UTF-8
#'
#' @import dplyr
#' @import sf
NULL
