#' Join bands of two identically shaped data cubes 
#' 
#' Create a proxy data cube, which joins the bands of two identically shaped data cubes. The resulting cube
#' will have bands from both input cubes.
#'
#' @param cube_list a list with two or more source data cubes
#' @param cube_names list or character vector with optional name prefixes for bands in the output data cube (see Details)
#' @return proxy data cube object
#' @examples 
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
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-05"),
#'                           srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube.b04 = select_bands(raster_cube(L8.col, v), c("B04"))
#' L8.cube.b05 = select_bands(raster_cube(L8.col, v), c("B05"))
#' join_bands(list(L8.cube.b04,L8.cube.b05))
#' \donttest{
#' plot(join_bands(list(L8.cube.b04,L8.cube.b05)))
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details 
#' The number of provided cube_names must match the number of provided input cubes.
#' If no cube_names are provided, bands of the output cube will adopt original names from the input cubes (without any prefix). If any two of the input bands have identical names,
#' prefixes default prefixes ("X1", "X2", ...) will be used.
#' 
#' @export 
join_bands <- function(cube_list, cube_names = NULL) {
  
  
  ncubes = length(cube_list)
  
  stopifnot(is.list(cube_list))
  stopifnot(length(cube_list) > 1)

  for (i in 1:ncubes) {
    stopifnot(is.cube(cube_list[[i]]))
  }

  if (is.null(cube_names)) {
    names = character()
    for (i in 1:ncubes) {
      names = c(names, names(cube_list[i])) 
    }
    if (anyDuplicated(names) > 0) {
      cube_names = paste0("X", 1:ncubes)
    }
    else {
      cube_names = character()
    }
  }
  else {
    if (is.list(cube_names)) {
      cube_names = as.character(unlist(cube_names))
    }
    stopifnot(length(cube_names) == ncubes)
    stopifnot(is.vector(cube_names))
    stopifnot(is.character(cube_names))
  }
  
  
  
  x = libgdalcubes_create_join_bands_cube(cube_list, cube_names)
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




