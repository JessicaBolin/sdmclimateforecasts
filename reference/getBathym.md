# Extract bathymetry

This function extracts bathymetry for lat/lons of interest directly from
ERDAPP.

## Usage

``` r
getBathym(points, desired.diameter = 0.7, func = "mean")
```

## Arguments

- points:

  Data frame. Input data frame with columns labeled `lon` and `lat`

- desired.diameter:

  Numeric. Defaults to 0.7˚. Results in a 7x7 (= 49) pixel box
  (0.7x0.7˚) around the central location.

- func:

  Character. Default "mean", can change to "sd".

## Value

A five column data frame with bathymetry extracted for lat/lons of
interest.

## Details

Can be slow for large datasets.

## Examples

``` r
points <- expand.grid(
  lon = seq(-132, -122, by = 2),
  lat = seq(32, 40, by = 2)
)

testbathym <- getBathym(
  points = points,
  desired.diameter = 0.7,
  func = "mean"
)
#> Registered S3 method overwritten by 'httr':
#>   method           from  
#>   print.cache_info hoardr
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%

head(testbathym)
#>    lon lat bathym_mean_0.7 bathym_count_0.7
#> 1 -132  32       -4550.447             1849
#> 2 -130  32       -4502.723             1849
#> 3 -128  32       -4279.393             1849
#> 4 -126  32       -4207.784             1849
#> 5 -124  32       -4198.538             1849
#> 6 -122  32       -4043.458             1849
```
