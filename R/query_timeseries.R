
#' Query data cube timeseries at irregular spatial points
#' 
#' This function will overlay provided spatial points with a data cube and return time series of all bands
#' of the cells for each query point, as a list of data.frame (one data frame per band) where rows correspond to points and columns represent time.
#' If needed, point coordinates are automatically transformed to the SRS of the data cube.
#'
#' @param x source data cube
#' @param px vector of x coordinates
#' @param py vector of y coordinates
#' @param srs spatial reference system string identifer (as GDAL understands) 
#' @return a list of data.frames (one per band / variable) with one row per point and one column per data cube time slice
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
#' x = seq(from = 388941.2, to = 766552.4, length.out = 10)
#' y = seq(from = 4345299, to = 4744931, length.out = 10)
#' 
#' query_timeseries(L8.rgb, x, y, srs(L8.rgb))
#' 
#' @export
query_timeseries <- function(x, px, py, srs) {
  if (length(px) != length(py)) {
    stop("Expected identical length for point coordinates px and py.")
  }
  
  res = libgdalcubes_query_timeseries(x, px, py, srs)
  ts = lapply(res, function(z) {
    z = as.data.frame(z)
    colnames(z) = dimension_values(x)$t
    z
  })
  names(ts) = names(x)
  return(ts)
}