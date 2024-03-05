# gdalcubes 0.7.0 (2024-03-06)

* add `as.data.frame()` to easily convert data cubes to data frames
* add `predict.cube()` to predict pixel values based on models (lm, glm, caret, tidymodels, and similar)
* add `window_space()` to apply (focal) moving window kernels or aggregation functions
* `extract()` now combines extracted values with input geometries and attributes (if `merge = TRUE`)           
* add support for imagery with spatial reference from geolocation arrays (including curvilinear grids) 
* `stac_image_collection()` now accepts STACItemCollection objects directly and should be more robust
* Windows build uses pkg-config if available
* Improved error reporting for inaccessible imagery

# gdalcubes 0.6.4 (2023-04-14)

* add native quartile reducers in `reduce_time()`
* fix r-devel UCRT win build on CRAN
* fix crashes on Windows UCRT due to unusable std::regex()
* fix parallel data cube processing when nonstandard external package locations are used
* `stack_cube()` now ignores files if not accessible / invalid instead of stopping all computations
* The codebase has been reorganized R package is now maintained under https://github.com/appelmar/gdalcubes, whereas the C++ repo will be archived. 

# gdalcubes 0.6.3 (2023-01-19)

* fix gcc-13 compiler errors on CRAN
* add datetime interval support in STAC collections
* add support of new windows toolchain using Makevars.ucrt


# gdalcubes 0.6.2 (2022-10-09)

* fix clang-15 compiler warnings on CRAN
* new operation `aggregate_space()` to reduce spatial resolution of data cubes
* improved / faster implementation of `plot()` 
* handle WKT strings as spatial reference systems in STAC responses
* handle special characters in variable / band names


# gdalcubes 0.6.1 (2022-03-22)

* fix gcc-12 builds on CRAN
* fix automatic reprojection in `extract_geom()`
* update GDAL on Windows 

# gdalcubes 0.6.0 (2022-03-07)

* major stability improvements:
  * fix unexpected stack overflows due to to GDAL error handler from `sf` calling `Rf_warning()`
  * if GDALOpen() fails to read an image, it will now be simply ignored but not stop processing the current chunk
  * improved handling and checks for empty chunks in data cube operations 
  * parallel processing now uses worker processes instead of threads
* new `extract_geom()` function to extract data cube values from spatial or spatiotemporal features and to compute summary statistics
* remove functions `query_points()`, `query_timeseries()`, and `zonal_statistics()` in favor of `extract_geom()`
* fix `filter_geom()` issues with larger polygons 
* fix `filter_geom()` error while checking if polygon is within data cube
* use WKT strings or authority codes in image collections instead of proj4 strings
* chunk sizes can now be set as a global package option either as constant sizes or as a function of data cube size
* default chunk sizes consider the number of parallel worker processes 
* `animate()` now can produce mp4 and GIF animations
* `animate()` works for larger image sequences using the `av` or `gifski` packages
* remove dependency on `RcppProgress`


# gdalcubes 0.5.1 (2021-02-12)

* fix CRAN vignette issue on Mac due to data download failures
* fix `image_mask()` function for minimum and maximum values 

# gdalcubes 0.5.0 (2021-10-27)

* new operation `aggregate_time()` to reduce temporal resolution of data cubes
* new `stack_cube()` function to build data cubes from aligned images without image collection creation
* new operations `slice_time()` and `slice_space()` to extract single time series or slices
* new `crop()` function can be used to crop a data cube by space and/or time
* single bands of data cubes can be selected using the `$` operator
* fix datetime parser to support strings with fractional seconds
* fix CRAN issues due to obsolete autoconf warnings
* the`[]` operator can now be used for flexible cropping, slicing, and band selection on data cubes


# gdalcubes 0.4.1 (2021-07-29)

* fix build issues on MacOS
* fixes for Windows build including ucrt support


# gdalcubes 0.4.0 (2021-07-08)

* new operator `ncdf_cube()` to read data cubes from (intermediate) results
* new operator `rename_bands()` to change band names
* image collection creation from STAC API queries with `stac_image_collection()`
* progress bar can now be disabled with `gdalcubes_options()`
* removed `gdalcubes_set_threads()` in favor of `gdalcubes_options()` 
* removed `gdalcubes_debug_output()` in favor of `gdalcubes_options()` 
* removed `gdalcubes_set_ncdf_compression()` in favor of `gdalcubes_options()` 
* removed `gdalcubes_use_cache()` in favor of `gdalcubes_options()` 
* removed image collection operations `translate_COG()` and `translate_gtiff()`
* fix installation issues on MacOS and GCC11 warnings


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