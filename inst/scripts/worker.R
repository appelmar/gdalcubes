#!/usr/bin/env r

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Missing input argument(s)")
} 

json = args[1]
#json = ".test_worker_description.json" # for testing

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("package jsonlite required; please install first")
}
j = jsonlite::fromJSON(json)


library(gdalcubes)
# set gdalcubes options
do.call(gdalcubes_options,args = j$gdalcubes_options)

# set GDAL config options
gdal_options = j$gdal_options
if (!is.null(gdal_options)) {
  for (key in names(gdal_options)) {
    gdalcubes_set_gdal_config(key, gdal_options[key])
  }
}

gdalcubes:::gc_exec_worker(j$cube, j$worker_id, j$worker_count, j$workdir, gdalcubes_options()$ncdf_compression_level)
  