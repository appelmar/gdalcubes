
#' Query data cube values at irregular spatiotemporal points
#' 
#' This function will overlay provided spatiotemporal points with a data cube and return all band values
#' of the cells for each query point, as a data.frame where rows correspond to points and columns represent bands.
#' If needed, point coordinates are automatically transformed to the SRS of the data cube.
#'
#' @param x source data cube
#' @param px vector of x coordinates
#' @param py vector of y coordinates
#' @param pt vector of date/ time coordinates 
#' @param srs spatial reference system string identifer (as GDAL understands) 
#' @return a data.frame with one row per point and one column per data cube band or variable
#' @details 
#' 
#' Date and time of the query points can be provided as vector of class character, Date, or POSIXct.
#' 
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
#' t = seq(as.Date("2018-01-01"), as.Date("2018-12-02"), length.out = 10 )
#' 
#' query_points(L8.rgb, x, y, t, srs(L8.rgb))
#' 
#' @export
query_points <- function(x, px, py, pt, srs) {
  
  
  if (length(px) != length(py) || length(py) != length(pt)) {
    stop("Expected identical length for point coordinates px, py, and pt.")
  }
  
  if (!is.character(pt)) {
    pt = format(pt)
  }
  
  df = as.data.frame(libgdalcubes_query_points(x, px, py, pt, srs))
  colnames(df) <- names(x)
  return(df)
  
}