
#' Convert a data cube to an in-memory R array
#' 
#' @param x data cube
#' @return Four dimensional array with dimensions band, t, y, x
#' @examples 
#' \donttest{
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-05"),
#'               srs="EPSG:32618", nx = 100, ny=100, dt="P1M")
#' as_array(select_bands(raster_cube(L8.col, v), c("B04", "B05")))
#' }
#' @note Depending on the data cube size, this function may require substantial amounts of main memory, i.e.
#' it makes sense for small data cubes only.
#' @export
as_array <- function(x) {
  
  stopifnot(is.cube(x))
  size = c(nbands(x), size(x))
  
  
  fn = tempfile(fileext = ".nc")
  libgdalcubes_eval_cube(x, fn, .pkgenv$compression_level)
  
  
  f <- ncdf4::nc_open(fn)
  # derive name of variables but ignore non three-dimensional variables (e.g. crs)
  vars <- names(which(sapply(f$var, function(v) {
    if (v$ndims == 3)
      return(v$name)
    return("")
  }) != ""))
  
  out = array(NA, dim=size)
  for (b in 1:length(vars)) {
    z = aperm(ncdf4::ncvar_get(f, vars[b], collapse_degen=FALSE), c(3,2,1))  # does it drop dims?
    dim(z) <- size[2:4]
    out[b,,,] = z
  }
  ncdf4::nc_close(f)
  
  dv <- dimension_values(x)
  dimnames(out) <- list(bands=vars, t=dv$t, y=dv$y, x=dv$x)
  
  return(out)
}



