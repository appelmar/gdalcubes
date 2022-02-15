#' Create a dummy data cube with a fill value
#' 
#' Create a data cube with a constant fill value for one or more bands from a data cube view. Use this cube for testing.
#' 
#' @param view a data cube view defining the shape (spatiotemporal extent, resolution, and spatial reference)
#' @param chunking length-3 vector or a function returning a vector of length 3, defining the size of data cube chunks in the order time, y, x
#' @param fill fill value
#' @param nbands number of bands 
#' @details 
#' If chunking is provided as a function, it must accept exactly three arguments for the total size of the cube in t, y, and x axes (in this order).
#' 
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
#' @noRd
.raster_cube_dummy <- function(view, nbands=1, fill=1, chunking=.pkgenv$default_chunksize) {
  if (missing(view)) {
    stop("missing view argument")
  }
  if (is.function(chunking)) {
    chunking = chunking(view$time$nt, view$space$ny, view$space$nx)
  }
  x = gc_create_dummy_cube(view, nbands, fill, chunking)
  class(x) <- c("dummy_cube", "cube", "xptr")
  return(x)
}


#' Create a completely empty dummy data cube 
#' 
#' This cube can be used to test whether data cube operations 
#' correctly deal with empty chunks. 
#' 
#' @param view a data cube view defining the shape (spatiotemporal extent, resolution, and spatial reference)
#' @param chunking length-3 vector or a function returning a vector of length 3, defining the size of data cube chunks in the order time, y, x
#' @param nbands number of bands 
#' @details 
#' If chunking is provided as a function, it must accept exactly three arguments for the total size of the cube in t, y, and x axes (in this order).
#' 
#' @return a proxy data cube object
#' @examples 
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'                          bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'              srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = .raster_cube_empty(v, 1)
#' L8.cube
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @noRd
.raster_cube_empty <- function(view, nbands=1, chunking=.pkgenv$default_chunksize) {
  if (missing(view)) {
    stop("missing view argument")
  }
  if (is.function(chunking)) {
    chunking = chunking(view$time$nt, view$space$ny, view$space$nx)
  }
  x = gc_create_empty_cube(view, nbands, chunking)
  class(x) <- c("empty_cube", "cube", "xptr")
  return(x)
}







