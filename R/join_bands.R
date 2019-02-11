#' Join bands of two identically shaped data cubes 
#' 
#' Create a proxy data cube, which joins the bands of two identically shaped data cubes. The resulting cube
#' will have bands from both input cubes.
#'
#' @param A first source data cube
#' @param B second source data cube
#' @return proxy data cube object
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"), 
#'                         ".TIF", recursive = TRUE, full.names = TRUE) 
#' 
#' v = cube_view(extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'                           srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.col = create_image_collection(L8_files, "L8_L1TP")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube.b04 = select_bands(raster_cube(L8.col, v), c("B04"))
#' L8.cube.b05 = select_bands(raster_cube(L8.col, v), c("B05"))
#' join_bands(L8.cube.b04,L8.cube.b05)
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export 
join_bands <- function(A, B) {
  stopifnot(is.cube(A))
  stopifnot(is.cube(B))
  
  x = libgdalcubes_create_join_bands_cube(A, B)
  class(x) <- c("join_bands_cube", "cube", "xptr")
  return(x)
}



is.join_bands_cube  <- function(obj) {
  if(!("join_bands_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




