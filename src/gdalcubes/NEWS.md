
# 0.3.0

* New data cube types `select_time_cube`, `filter_geom_cube`, `stream_apply_time_cube`
* New data cube extraction functions `query_timeseries` and `zonal_statistics`
* New batch functions to translate complete imagery of image collections to GeoTIFF / COG 
* Revised data model to support irregular / labeled time dimension
* Replaced external JSON library with `json11`
* New gdalwarp client to avoid multithreading issues with GDAL3
* fix CRS metadata in produced netCDF files
* Optional global SRS definition in collection formats
* `join_bands_cube` may now combine more than two data cubes

# 0.2.3

* fixed clang compiler warnings
* fixed MODIS collection formats
* new collection formats MxD14A2 and MxD13A2


# 0.2.2

* support for GDAL subdatasets in collection formats
* add `query_points` operation to query data cube values at irregular spatiotemporal points


# 0.2.1

* fixed compiler warnings
* fix mean aggregation


# 0.2.0

* reimplementation of `image_collection_cube::read_chunk()`, supporting image masks and custom gdalwarp arguments
* renamed `filter_predicate` -> `filter_pixel`
* added gdalcubes namespace
* removed reduce and stream commands from command line interface
* added `gdalcubes exec` command to command line interface
* added experimental `gdalcubes exec` and `gdalcubes translate_cog` for batch processing of image collections
* added `stream_reduce_time operator`, applying external processes on time series
* added `keep_bands` option to `apply_pixel` operation, keeping all bands from the input data cube


# 0.1.1

* removed exprtk library
* output NetCDF files now contain bounds variables
* readthedocs theme for documentation
* gdalwarp now receives correct extent (in full double precision)

# 0.1.0

* First Release