# gdalcubes 0.2.3 (2019-10-21)

* fixed clang compiler warnings 
* fixed MODIS collection formats
* new collection formats MxD14A2 and MxD13A2

# gdalcubes 0.2.2 (2019-10-15)

* support for GDAL subdatasets in collection format
* MODIS collection formats now use subdatasets automatically
* fixed configure.ac for R-devel
* add `query_points()` to query data cube values at irregular spatiotemporal points


# gdalcubes 0.2.1 (2019-08-21)

* new collection format for PlanetScope data
* fixed R CMD check warnings on CRAN (caused by compiler warning -Wdeprecated-declarations)
* fixed mean aggregation


# gdalcubes 0.2.0 (2019-08-07)

## New Features
* add `animate()` function to create data cube time series animations
* apply mask bands on pixel values during the construction of the data cube, see `?image_mask`
* add `write_tif()` to export data cubes as (possibly cloud-optimized) GeoTIFF files (one per time slice)
* export of data cubes with `write_tif()` and `write_ncdf()` supports packing data values to smaller integer types  
* processing cubes is interruptible, though it can still take time to let all threads finish their current chunk
* add `as_array()` function to convert a data cube to a native in-memory R array
* new operator `fill_time()` fills NA pixels of data cubes based on time series interpolation
* changed image collection database schema, existing collections must be recreated
* new global configuration function `gdalcubes_options()` as a replacement to `gdalcubes_set_threads()` etc.
* new function `add_images()` adds images to an existing image collection

## Minor improvements
* rename `filter_predicate()` -> `filter_pixel()`
* collection format Sentinel2_L2A now includes WVP, AOT, and SCL bands 
* consistent output for printing data cube views and data cubes
* new collection format for Sentinel-2 data on Theia (credits to Xavier Laviron)
* new collection format for MODIS MxD13Q1 vegetation index data
* add `write_json_descr`argument to `write_ncdf()`
* new argument `with_VRT` in `write_ncdf()` to write GDAL VRT datasets for data cube time slices
* collection formats can now overwrite scale, offset, and unit for bands
* `write_ncdf()` can produce netCDF files without bounds variables if desired
* `write_ncdf()` and `write_tif()` return created files as character vectors.

## Bug fixes
* fix windows source compilation on CRAN
* bands of multiband files are now read in correct order
* fix package build with PROJ 6.1 (credits to Roger Bivand)




# gdalcubes 0.1.0 (2019-05-15)

* First release