#' Extract bathymetry
#'
#'This function extracts bathymetry for lat/lons of interest directly from ERDAPP.
#'
#' @param points Data frame. Input data frame with columns labeled `lon` and `lat`
#' @param desired.diameter Numeric. Defaults to 0.7˚. Results in a 7x7 (= 49) pixel box (0.7x0.7˚) around the central location.
#' @param func Character. Default "mean", can change to "sd".
#'
#' @return A five column data frame with bathymetry extracted for lat/lons of interest.
#'
#' @details
#' Can be slow for large datasets.
#'
#' @export
#'
#' @importFrom rerddap
#' @importFrom rerddapXtracto
#'
#' @examples
#' points <- expand.grid(
#'   lon = seq(-132, -122, by = 2),
#'   lat = seq(32, 40, by = 2)
#' )
#'
#' testbathym <- getBathym(
#'   points = points,
#'   desired.diameter = 0.7,
#'   func = "mean"
#' )
#'
#' head(testbathym)


getBathym <- function(points, desired.diameter = 0.7, func = "mean") {

  fishlat <- points$lat
  # See what longitude is called, and whether it's in degrees east or west
  # I often use "lon360" to specify longitude in degrees east: if that's the case, just grab that
  if("lon360" %in% colnames(points)) {
    fishlon360 <- points$lon360
  } else if("lon" %in% colnames(points)) { # Or maybe there's a column called "lon", and we don't know how it's measured
    fishlon360 <- points$lon
  }
  # Convert to degrees east if needed
  fishlon360 <- ifelse(fishlon360 < 0, fishlon360 + 360, fishlon360)

  # Get dataset info from ERDDAP
  bathInfo <- rerddap::info('etopo360', url = 'http://coastwatch.pfeg.noaa.gov/erddap/')
  # Extract bathymetry at points
  bathym <- rerddapXtracto::rxtracto(bathInfo, parameter = 'altitude', xcoord = fishlon360, ycoord = fishlat,
                     xlen = desired.diameter, ylen = desired.diameter, progress_bar = TRUE)
  bathym <- data.frame(bathym)

  # Create cols
  points$out <- NA

  # Output values
  if(func == "mean") {
    points$out <- bathym$mean.altitude
    colnames(points)[ncol(points)] <- paste0("bathym_mean", "_", desired.diameter)
  } else if (func == "sd") {
    points$out <- bathym$stdev.altitude
    colnames(points)[ncol(points)] <- paste0("bathym_sd", "_", desired.diameter)
  }
  points$count <- bathym$n
  colnames(points)[ncol(points)] <- paste0("bathym_count", "_", desired.diameter)
  # Return
  return(points)
}
