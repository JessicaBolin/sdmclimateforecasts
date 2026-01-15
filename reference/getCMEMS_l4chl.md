# Extract CMEMS level 4 daily chl interpolated product

Requires global 4km netcdfs downloaded from the [CMEMS L4 global ocean
color
product](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description)
and stored locally.

## Usage

``` r
getCMEMS_l4chl(points, desired.diameter = 0.7, func = c("mean", "sd"), nc.path)
```

## Arguments

- points:

  Data frame. Requires an input data frame with columns labeled "lon"
  (degrees west), "lat", and "date" (formatted as dates).

- desired.diameter:

  Numeric. Default = 0.7. Results in a 0.7 x 0.7 degree box around the
  central location

- func:

  Character. "mean" or "sd"

- nc.path:

  Character. File path to where daily netCDF's are stored

## Value

A five column data frame with chla extracted for lat/lons/dates of
interest.

## Details

Unlike ROMS, CMEMS chl resolution is defined in km, not degrees. So
pixel sizes will vary slightly with latitude. To keep workflow
consistent with existing FRD/ESD workflows, we assume that pixel sizes
are 0.0416667 X 0.0416667. This may not be accurate at high latitudes!

## Examples

``` r
points <- data.frame("lon" = c(-130, -125, -124, -130, -125, -124, -130, -125, -124, -180),
"lat" = c(45, 34, 32, 40, 45, 34, 32, 40, 36, 41),
"date" = as.Date("2011-04-04"))
desired.diameter <- 0.7 # Note! This means a 0.7x0.7 degree box
func <- "mean" # mean or sd
nc.path <- "/Users/admin/Downloads"
#getCMEMS_l4chl(points, desired.diameter = 0.7, "mean", nc.path)
```
