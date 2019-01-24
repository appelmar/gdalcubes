
#' Apply arithmetic expressions over all pixels of a data cube
#' 
#' Create a proxy data cube, which applies arithmetics expressions over all pixels of a data cube. Expressions may access band values by their name.
#'
#' @param cube Source data cube
#' @param expr character vector with one or more arithmetic expressions (see Details)
#' @param names optional character vector with the same length as expr to specify band names for the output cube
#' @return A proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/ArashPartow/exprtk}{exprtk library} to evaluate expressions in C++, you can look at the library examples 
#' to see what kind of expressions you can execute.
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"), 
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = gcbs_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = gcbs_create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = gcbs_cube(L8.col, v) 
#'  L8.cube = gcbs_select_bands(L8.cube, c("B04", "B05")) 
#'  L8.ndvi = gcbs_apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#'  L8.ndvi
#'  L8.ndvi.median =  gcbs_reduce_time(L8.ndvi, "median(NDVI)") 
#'  plot(L8.ndvi.median, key.pos=1, zlim=c(0,1))
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
gcbs_apply_pixel <- function(cube, expr, names=NULL) {
  stopifnot(is.gcbs_cube(cube))
  
  if (is.null(names)) {
    names <- paste("band", 1:length(expr), sep="")
  }
  
  x = libgdalcubes_create_apply_pixel_cube(cube, expr, names)
  class(x) <- c("gcbs_apply_pixel_cube", "gcbs_cube", "xptr")
  return(x)
}



is.gcbs_apply_pixel_cube  <- function(obj) {
  if(!("gcbs_apply_pixel_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




