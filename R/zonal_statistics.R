
#' Query summary statistics of data cube values over polygons
#' 
#' This function will overlay spatial polygons with a data cube and compute summary statistics
#' of pixels within polygons over all time slices. 
#'
#' @param x source data cube
#' @param geom Either an sf object, or a path to an OGR dataset (Shapefile, GeoPackage, or similar) with input polygon geometries 
#' @param expr character vector of summary statistics expressions, describing pairs of aggregation functions and data cube bands (e.g. "mean(band1)")  
#' @param out_dir output directory where resulting files will be written to
#' @param prefix prefix of resulting filenames; will be followed by date/time
#' @param ogr_layer If the input OGR dataset has multiple layers, a layer can be chosen by name
#' @return vector of paths to produced GeoPackage files 
#' @details 
#' 
#' Input geometries must be stored as OGR dataset (e.g. shapefile, geopackage, or similar). 
#' 
#' The function creates one geopackage per time slice. Each geopackage contains all polygons as geometries and all computed 
#' statistics as attributes. 
#' 
#' Available summary statistics currently include "min", "max", "mean", "median", "count", "sum", and "prod". 
#' 
#' 
#' @note Currently, the spatial reference systems of the data cube and the feautures must be identical.
#' 
#' 
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01-01", t1="2018-12-02"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P14D")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' 
#' #TODO
#' 
#' @export
zonal_statistics <- function(x, geom, expr, out_dir = tempdir(), prefix=basename(tempfile()), ogr_layer=NULL) {

  
  stopifnot(is.cube(x))
  
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist")
  }
  
  if ("sf" %in% class(geom)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("sf package not found, please install first") 
    geom_file = tempfile(fileext = ".gpkg")
    st_write(st_geometry(geom), geom_file, quiet = TRUE)
    geom = geom_file
  }
  else if ("sfc" %in% class(geom)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("sf package not found, please install first") 
    geom_file = tempfile(fileext = ".gpkg")
    st_write(geom, geom_file, quiet = TRUE)
    geom = geom_file
    
  }
  if (!is.character(geom) || length(geom) > 1) {
    stop("geom must be either a character vector with length 1, or an sf object")
  }
  
  
  if (is.null(ogr_layer)) {
    ogr_layer = ""
  }
   
  agg_funcs = gsub("\\(.*\\)", "", expr)
  agg_bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(agg_funcs) == length(agg_bands))
  
  return(libgdalcubes_zonal_statistics(x, geom, agg_funcs, agg_bands, out_dir, prefix, ogr_layer))
  
}