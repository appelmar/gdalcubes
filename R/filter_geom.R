
#' Filter data cube pixels by a polygon
#' 
#' Create a proxy data cube, which filters pixels by a spatial (multi)polygon For all pixels whose center is within the polygon, the original
#'
#' @param cube source data cube
#' @param geom either a WKT string, or an sfc or sfg object (sf package)
#' @param srs string identifier of the polygon's coordinate reference system understandable for GDAL
#' @return a proxy data cube object
#' @details The resulting data cube will not be cropped but pixels outside of the 
#' polygon will be set to NAN. 
#' 
#' If \code{geom} is provided as an sfc object with length > 1, geometries will
#' be combined with \code{sf::st_combine()} before.
#' 
#' The geometry is automatically transformed to the data cube's spatial reference
#' system if needed.
#'  
#' @examples
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI")
#' WKT = gsub(pattern='\\n',replacement="",x = 
#'   "Polygon ((-74.3541 40.9254, 
#'              -73.9813 41.2467, 
#'              -73.9997 41.4400, 
#'              -74.5362 41.1795, 
#'              -74.6286 40.9137, 
#'              -74.3541 40.9254))")
#' L8.ndvi.filtered = filter_geom(L8.ndvi, WKT, "EPSG:4326")
#' L8.ndvi.filtered
#' \donttest{
#' plot(L8.ndvi.filtered)
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
filter_geom <- function(cube, geom, srs = NULL) {
  stopifnot(is.cube(cube))
  
  if ("sfc" %in% class(geom) || "sfg" %in% class(geom)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("package sf required; please install first")
    }
    if (is.null(srs)) {
      temp = sf::st_crs(geom)
      if(!is.null(temp$wkt)) {
        srs = temp$wkt
      }
      else if(!is.null(temp$proj4string)) {
        srs = temp$proj4string
      }
      else {
        stop ("cannot fetch coordinate reference system of geometry")
      }
    }
    if (length(geom) > 1) {
      geom = sf::st_combine(geom)
    }
    geom = sf::st_as_text(geom)
  }
  else if (is.character(geom)) {
    if (is.null(srs)) {
      stop("srs argument required for WKT input")
    }
  }
 
  x = libgdalcubes_create_filter_geom_cube(cube, geom, srs)
  class(x) <- c("filter_geom_cube", "cube", "xptr")
  return(x)
}



is.filter_geom_cube  <- function(obj) {
  if(!("filter_geom_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




