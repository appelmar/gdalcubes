# gdalcubes 0.1.9999

* rename `filter_predicate()` -> `filter_pixel()`
* new operator `fill_time()` to fill NA pixels of data cubes based on time series interpolation
* add masking based on pixel band values while reading images, see `?image_mask`
* collection format Sentinel2_L2A now includes WVP, AOT, and SCL bands 
* add `write_json_descr`argument to `write_ncdf()`
* `write_ncdf`, `as_stars()`, and `plot()` are now  interruptible, though it might take some time to let all threads finish their current chunk


# gdalcubes 0.1.0 (2019-05-15)

* First release