#' Select bands of a data cube
#' 
#' Create a proxy data cube, which selects specific bands of a data cube. The resulting cube
#' will drop any other bands.
#'
#' @param cube Input data cube
#' @param bands character vector with band names
#' @return A proxy data cube object
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = cube_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = data_cube(L8.col, v) 
#'  L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#'  L8.rgb
#'  L8.rgb.median = reduce_time(L8.rgb, "median(B02)", "median(B03)", "median(B04)")  
#'  plot(L8.rgb.median, rgb=3:1, zlim=c(4000,12000))
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @note For performance reasons, select_bands should always be called directly on a image collection cube created with cube and 
#' drop all unneded bands.This allows to reduce GDAL RasterIO and warp operations.
#' @export
select_bands <- function(cube, bands) {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_select_bands_cube(cube, bands)
  class(x) <- c("select_bands_cube", "cube", "xptr")
  return(x)
}




is.select_bands_cube  <- function(obj) {
  if(!("select_bands_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




