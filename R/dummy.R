
#' Create a dummy data cube with a fill value
#' 
#' Create a data cube with a constant fill value for one or more bands from a data cube view. Use this cube for testing.
#' 
#' @param view a data cube view defining the shape (spatiotemporal extent, resolution, and spatial reference)
#' @param chunking vector of length 3 defining the size of data cube chunks in the order time, y, x.
#' @param fill fill value
#' @param nbands number of bands 
#' @return a proxy data cube object
#' @examples 
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'                          bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'              srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube_dummy(v, 1, 2.345)
#' \donttest{
#' plot(L8.cube, zlim=c(0,4))
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
raster_cube_dummy <- function(view, nbands=1, fill=1, chunking=c(16, 256, 256)) {
  x = libgdalcubes_create_dummy_cube(view, nbands, fill, chunking)
  class(x) <- c("dummy_cube", "cube", "xptr")
  return(x)
}

is.dummy_cube  <- function(obj) {
  if(!("dummy_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




