# changelog

> 2026-01-14: package creation and test with `getBathym`

- amended DESCRIPTION to include authors
- defined package dependencies for getBathym function
  (`use_package("__")`) - auto amends DESCRIPTION
- created function with `use_r("getBathym")` - inserted roxygen skeleton
  to build help file (insert cursor in script, then code -\> insert
  roxygen skeleton - Created `.rd` documentation with
  `devtools::document()`
- `usethis::use_mit_license()` - Added MIT licence
- Added some dummy vignettes
  - `use_vignette("changelog.Rmd")`
- Added ‘Getting started’ page
- Updated `_pkgdown.yml`
- Added home index landing page `README.md`
