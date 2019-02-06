
#' Filter data cube pixels by a user-defined predicate on band values
#' 
#' Create a proxy data cube, which evaluates a predicate over all pixels of a data cube. For all pixels which fulfill the predicate, the original
#' band values are returned. Other pixels are simply filled with NANs. The predicate may access band values by their name.
#'
#' @param cube Source data cube
#' @param pred predicate to be evaluated over all pixels
#' @return A proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/ArashPartow/exprtk}{exprtk library} to evaluate expressions in C++, you can look at the library examples 
#' to see what kind of expressions you can execute.
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'                bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'                srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#'  L8.col = create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = data_cube(L8.col, v) 
#'  L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#'  L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#'  L8.ndvi.filtered = filter_predicate(L8.ndvi, "NDVI > 0.5") 
#'  L8.ndvi.filtered
#'  L8.ndvi.count =  reduce_time(L8.ndvi, "count(NDVI)") 
#'  plot(L8.ndvi.count, key.pos=1)
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
filter_predicate <- function(cube, pred) {
  stopifnot(is.cube(cube))

  x = libgdalcubes_create_filter_predicate_cube(cube, pred)
  class(x) <- c("filter_predicate_cube", "cube", "xptr")
  return(x)
}



is.filter_predicate_cube  <- function(obj) {
  if(!("filter_predicate_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




