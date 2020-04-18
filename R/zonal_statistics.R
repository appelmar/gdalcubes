
#' Query summary statistics of data cube values over polygons
#' 
#' This function will overlay spatial polygons with a data cube and compute summary statistics
#' of pixel values within the polygons over all time slices. 
#'
#' @param x source data cube
#' @param geom Either an sf object, or a path to an OGR dataset (Shapefile, GeoPackage, or similar) with input (multi)polygon geometries 
#' @param expr character vector of summary statistics expressions, describing pairs of aggregation functions and data cube bands (e.g. "mean(band1)")  
#' @param out_path path to where resulting GeoPackage will be written to
#' @param overwrite logical; overwrite \code{out_path} if file already exists, defaults to FALSE
#' @param ogr_layer If the input OGR dataset has multiple layers, a layer can be chosen by name
#' @return character length-one vector containing the path to the resulting GeoPackage file (see Details)
#' @details 
#' 
#' The function creates a single GeoPackage output file containing:
#' \itemize{
#'    \item A single layer "geom" containing the geometries (and feature identifiers) only.
#'    \item Attribute tables (layers without geometry) for each time slice of the data cube containing summary statistics as columns. 
#'    Corresponding layer names start with "attr_", followed by date and time.
#'    \item Virtual spatial views for each time slice, joining the geometries and attribute tables. 
#'    Corresponding layer names start with "map_", followed by date and time.
#' }
#' You will most-likely want to use the spatial view layers directly e.g. with the sf package.  
#' 
#' 
#' Available summary statistics currently include "min", "max", "mean", "median", "count", "sum", "prod", "var", and "sd". 
#' 
#' 
#' @note Currently, the spatial reference systems of the data cube and the features must be identical.
#' 
#' 
#' @examples 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"))
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(srs="EPSG:32618", dy=300, dx=300, dt="P1M", 
#'               aggregation = "median", resampling = "bilinear",
#'               extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931, 
#'                           t0="2018-01-01", t1="2018-12-31"))
#' L8.cube = raster_cube(L8.col, v) 
#' L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#' L8.ndvi
#' 
#' # toy example: overlay NDVI data with NYC districts
#' x = zonal_statistics(L8.ndvi, system.file("nycd.gpkg", package = "gdalcubes"),
#'                      expr = "median(NDVI)")
#' x
#' 
#' @export
zonal_statistics <- function(x, geom, expr, out_path = tempfile(fileext = ".gpkg"), overwrite = FALSE, ogr_layer = NULL) {

  
  stopifnot(is.cube(x))
  out_path = path.expand(out_path)
  if (file.exists(out_path) && !overwrite) {
    stop("Output file already exists; set overwrite = TRUE or choose a different output file")
  }
  
  if ("sf" %in% class(geom)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("sf package not found, please install first") 
    geom_file = tempfile(fileext = ".gpkg")
    sf::st_write(sf::st_geometry(geom), geom_file, quiet = TRUE)
    geom = geom_file
  }
  else if ("sfc" %in% class(geom)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("sf package not found, please install first") 
    geom_file = tempfile(fileext = ".gpkg")
    sf::st_write(geom, geom_file, quiet = TRUE)
    geom = geom_file
    
  }
  if (!is.character(geom) || length(geom) > 1) {
    stop("geom must be either a length-one character vector, or an sf object")
  }
  
  
  if (is.null(ogr_layer)) {
    ogr_layer = ""
  }
   
  agg_funcs = gsub("\\(.*\\)", "", expr)
  agg_bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(agg_funcs) == length(agg_bands))
  
  libgdalcubes_zonal_statistics(x, geom, agg_funcs, agg_bands, out_path, overwrite, ogr_layer)
  return(out_path)
}