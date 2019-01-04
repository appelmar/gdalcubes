
#' Set the number of threads for parallel data cube processing
#'
#' Data cubes can be processed in parallel where one thread processes one chunk at a time. Setting more threads
#' than the number of chunks of a cube thus has no effect and will not furhter reduce computation times.
#' @param n number of threads
#' @return nothing / invisible
#' @export
gcbs_set_threads <- function(n=1) {
  stopifnot(n >= 1)
  stopifnot(n%%1==0)
  libgdalcubes_set_threads(n)
  invisible()
}

#' Query gdalcubes version information
#'
#' @return List with gdalcubes library version information
#' @export
gcbs_version <- function() {
  return(libgdalcubes_version())
}

#' Enable debug output from the gdlcubes C++ library
#' @export
gcbs_enable_debug_output <- function() {
  libgdalcubes_debug_output(TRUE)
  invisible()
}

#' Disable debug output from the gdlcubes C++ library
#' @export
gcbs_disable_debug_output <- function() {
  libgdalcubes_debug_output(FALSE)
  invisible()
}
