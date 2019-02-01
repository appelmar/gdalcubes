#' Join bands of two identically shaped data cubes 
#' 
#' Create a proxy data cube, which joins the bands two identically shaped data cubes. The resulting cube
#' will have bands from both input cubes.
#'
#' @param A The first source data cube
#' @param B The second source data cube
#' @return A proxy data cube object
#' @examples 
#' \dontrun{join_bands(XXXX, YYYY)}
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export 
join_bands <- function(A, B) {
  stopifnot(is.cube(A))
  stopifnot(is.cube(B))
  
  x = libgdalcubes_create_join_bands_cube(A, B)
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




