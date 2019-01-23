
#' Filter data cube pixels by a user-defined predicate on band values
#' 
#' Create a proxy data cube, which evaluates a predicate over all pixels of a data cube. For all pixels which fulfill the predicate, the original
#' band values are returned. Other pixels are simply filled with NANs. The predicate may access band values by their name.
#'
#' @param cube Source data cube
#' @param pred predicate to be evaluated over all pixels
#' @return A proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/ArashPartow/exprtk}{exprtk library} to evaluate expressions in C++, you can look at the library examples 
#' to see what kind of expressions you can execute.
#' @examples 
#' \dontrun{
#' gcbs_filter_predicate(XXXX, expr=c("(B08-B04)/(B08+B04)"))
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
gcbs_filter_predicate <- function(cube, pred) {
  stopifnot(is.gcbs_cube(cube))

  x = libgdalcubes_create_filter_predicate_cube(cube, pred)
  class(x) <- c("gcbs_filter_predicate_cube", "gcbs_cube", "xptr")
  return(x)
}



is.gcbs_filter_predicate_cube  <- function(obj) {
  if(!("gcbs_filter_predicate_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




