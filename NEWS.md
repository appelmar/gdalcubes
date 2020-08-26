
# gdalcubes 0.3.1 (2020-08-25)

* make GEOS dependency optional
* remove CURL dependency from configure


# gdalcubes 0.3.0 (2020-08-05)

## New Features

* Compute summary statistics of data cubes over polygons with `zonal_statistics()` 
* Extracts time series at irregular spatial points with `query_timeseries()` 
* Time dimension may ow be irregular / labeled after selecting slices with the new `select_time()` function
* Filter pixels of a data cube by a spatial polygon with `filter_geom()`
* Apply an R function on time series without reduction using `apply_time()`
* Batch format conversion of images in a collection with `translate_cog()` and `translate_gtiff()`

## Minor improvements

* conversion to stars objects with `st_as_stars()`
* add support for image collections without collection format in `create_image_collection()`
* optional global SRS definition in collection formats
* default chunk size is now (t,y,x) = (1,256,256)
* remove `reduce()` function
* remove `cube` argument in `cube_view` function
* new collection format for daily 0.25° AVHRR Optimum Interpolation Sea Surface Temperature
* new collection formats for ESA CCI soil moisture products
* new collection format for daily precipitation observations from GPM / IMERG
* new collection format for MODIS MOD09GA (aqua and terra)
* add `na.color` argument in `plot.cube()`

## Bug fixes

* fix CRS metadata in produced netCDF files  
* fix multithreading locking issues with GDAL 3




# gdalcubes 0.2.5 (2020-05-17)

* fixed compiler warnings on CRAN
* temporarily removed `as_stars()`, will be added again in 0.3



# gdalcubes 0.2.4 (2020-02-02)

* fixed axis order issues with GDAL3 and PROJ6
* fixed compiler warnings with GDAL3



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