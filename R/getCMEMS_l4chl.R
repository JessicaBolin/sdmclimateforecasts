#' Extract CMEMS level 4 daily chl interpolated product
#'
#' Requires global 4km netcdfs downloaded from the [CMEMS L4 global ocean color product](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description) and stored locally.
#'
#' @param points Data frame.  Requires an input data frame with columns labeled "lon" (degrees west), "lat", and "date" (formatted as dates).
#' @param desired.diameter Numeric. Default = 0.7. Results in a 0.7 x 0.7 degree box around the central location
#' @param func Character. "mean" or "sd"
#' @param nc.path Character. File path to where daily netCDF's are stored
#'
#' @return A five column data frame with chla extracted for lat/lons/dates of interest.
#'
#' @details
#' Unlike ROMS, CMEMS chl resolution is defined in km, not degrees. So pixel sizes will vary slightly with latitude.
#' To keep workflow consistent with existing FRD/ESD workflows, we assume that pixel sizes are 0.0416667 X 0.0416667.
#' This may not be accurate at high latitudes!
#'
#' @export
#'
#' @examples
#' points <- data.frame("lon" = c(-130, -125, -124, -130, -125, -124, -130, -125, -124, -180),
#' "lat" = c(45, 34, 32, 40, 45, 34, 32, 40, 36, 41),
#' "date" = as.Date("2011-04-04"))
#' desired.diameter <- 0.7 # Note! This means a 0.7x0.7 degree box
#' func <- "mean" # mean or sd
#' nc.path <- "/Users/admin/Downloads"
#' #getCMEMS_l4chl(points, desired.diameter = 0.7, "mean", nc.path)



getCMEMS_l4chl <- function(points, desired.diameter = 0.7, func = c("mean", "sd"), nc.path) {

  # ---- Checks ----
  required_cols <- c("lon", "lat", "date")
  if (!all(required_cols %in% names(points))) {
    stop("Input points must have columns: lon, lat, date")
  }
  if (!inherits(points$date, "Date")) {
    stop("date column must be of class Date")
  }

  # ---- Resolution guard ----
  native_res <- 0.0416667
  desired.diameter <- max(desired.diameter, native_res)

  # ---- Longitude handling ----
  points$lon360 <- ifelse(points$lon < 0, points$lon + 360, points$lon)

  # ---- Output column ----
  stat <- match.arg(func, c("mean", "sd"))
  out_col <- paste0("chl_", stat, "_", desired.diameter)

  # remove existing column if present
  points[[out_col]] <- NULL
  points[[out_col]] <- NA_real_

  # ---- Date helpers ----
  fishyear  <- lubridate::year(points$date)
  fishmonth <- formatC(lubridate::month(points$date), width = 2, flag = "0")
  fishday   <- formatC(lubridate::day(points$date), width = 2, flag = "0")

  # ---- Loop ----
  for (i in seq_len(nrow(points))) {

    file1 <- file.path(
      nc.path,
      paste0(fishyear[i], fishmonth[i], fishday[i],
             "_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc")
    )

    file2 <- file.path(
      nc.path,
      paste0(fishyear[i], fishmonth[i], fishday[i],
             "_cmems_obs-oc_glo_bgc-plankton_myint_l4-gapfree-multi-4km_P1D.nc")
    )

    ncfile <- if (file.exists(file1)) file1 else if (file.exists(file2)) file2 else NA
    if (is.na(ncfile)) next

    nc <- ncdf4::nc_open(ncfile)
    on.exit(ncdf4::nc_close(nc), add = TRUE)

    lat <- ncdf4::ncvar_get(nc, "lat")
    lon <- ncdf4::ncvar_get(nc, "lon")
    lon360 <- ifelse(lon < 0, lon + 360, lon)

    if (is.na(points$lat[i]) || is.na(points$lon360[i])) next

    c0 <- which.min(abs(lon360 - points$lon360[i]))
    r0 <- which.min(abs(lat - points$lat[i]))

    pix_diam <- round(desired.diameter / native_res)
    pix_rad  <- floor(pix_diam / 2)

    c_low <- c0 - pix_rad
    c_up  <- c0 + pix_rad
    r_low <- r0 - pix_rad
    r_up  <- r0 + pix_rad

    ncol <- c_up - c_low + 1
    nrow <- r_up - r_low + 1

    if (c_low < 1) {
      v1 <- ncdf4::ncvar_get(nc, "CHL",
                             start = c(1, r_low, 1),
                             count = c(c_up, nrow, 1))
      v2 <- ncdf4::ncvar_get(nc, "CHL",
                             start = c(length(lon) + c_low, r_low, 1),
                             count = c(abs(c_low) + 1, nrow, 1))
      vals <- rbind(v1, v2)
    } else {
      vals <- ncdf4::ncvar_get(nc, "CHL",
                               start = c(c_low, r_low, 1),
                               count = c(ncol, nrow, 1))
    }

    if (all(is.na(vals))) next

    points[i, out_col] <-
      if (stat == "mean") mean(vals, na.rm = TRUE) else sd(vals, na.rm = TRUE)

    if (i %% 100 == 0) message(i, " points complete")
  }

  points
}

