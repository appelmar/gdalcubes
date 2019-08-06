# gdalcubes 0.2.0 (2019-08-06)

## New Features
* add `animate()` function to create data cube time series animations
* add masking based on pixel band values while reading images, see `?image_mask`
* add `as_array()` function to convert a data cube to a native in-memory R array
* add `write_tif()` to export data cube time slices as (possibly cloud-optimized) GeoTIFF files
* export of data cubes with `write_tif()` and `write_ncdf()` supports packing data values to smaller integer types  
* processing cubes is interruptible, though it might take some time to let all threads finish their current chunk
* new operator `fill_time()` fills NA pixels of data cubes based on time series interpolation
* changed image collection database schema, existing collections must be recreated
* new global configuration function `gdalcubes_options()` as a replacement to `gdalcubes_set_threads()` etc.

## Minor improvements
* rename `filter_predicate()` -> `filter_pixel()`
* collection format Sentinel2_L2A now includes WVP, AOT, and SCL bands 
* consistent output for printing data cube views and data cubes
* new collection format for Sentinel-2 data on Theia
* add `write_json_descr`argument to `write_ncdf()`
* new argument `with_VRT` in `write_ncdf()` to write GDAL VRT datasets for data cube time slices
* collection formats can now overwrite scale, offset, and unit for bands
* `write_ncdf()` can produce netCDF files without bounds variables if desired

## Bug fixes
* fix windows source compilation on CRAN
* bands of multiband files are now read in correct order






# gdalcubes 0.1.0 (2019-05-15)

* First release