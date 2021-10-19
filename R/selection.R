#' Select a data cube band by name
#' 
#' Select a data cube band by name
#' @name gdalcubes_selection
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P3M", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' L8.red = L8.cube$B04
#' 
#' \donttest{
#' plot(L8.red)
#' }
#' @param x source data cube
#' @param name character; name of selected band
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
"$.cube" = function(x, name) {
  stopifnot(is.cube(x))
  return(select_bands(x, name))
}
  