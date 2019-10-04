
#' Query data cube values at irregular spatiotemporal points
#' 
#' The function will overlay provided spatiotemporal points with the data cube and return all band values
#' of the cells per point, as a data.frame where rows correspoind to points and columns represent bands.
#' If needed, point coordinates are automativally converted to the SRS of the data cube.
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
#' TODO: what das GDAL unterstand?
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
#' TODO
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