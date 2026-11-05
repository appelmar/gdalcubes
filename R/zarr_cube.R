
#' Read a data cube from an existing Zarr dataset
#' 
#' Create a proxy data cube from a Zarr dataset that has been created using \code{\link{write_zarr}}.
#' This function does not read cubes from arbitrary Zarr datasets and can be used e.g., to load intermediate results and/or 
#' plotting existing cubes on disk without reprocessing the data cube creation, e.g. from image collections.
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' 
#' \donttest{
#' f = write_zarr(select_bands(raster_cube(L8.col, v), c("B02", "B03", "B04")))
#' zarr_cube(f)
#' }                           
#' 
#' @param path path to an existing Zarr dataset
#' @param chunking custom chunk sizes to read form the netCDF file; defaults to using chunk sizes from the netCDF file
#' @param auto_unpack logical; automatically apply offset and scale when reading data values
#' @return a proxy data cube object
#' @note This function returns a proxy object, i.e., it will delay any computations besides deriving the shape of the result.
#' @export
zarr_cube <- function(path, chunking = NULL, auto_unpack = TRUE) {
  stopifnot(file.exists(path))
  
  if (is.null(chunking)) {
    chunking = integer(0)
  }
  x = gc_create_zarr_cube(path, chunking, auto_unpack)
  class(x) <- c("zarr_cube", "cube", "xptr")
  return(x)
}



is.zarr_cube  <- function(obj) {
  if(!("zarr_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




