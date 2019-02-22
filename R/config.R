
#' Set the number of threads for parallel data cube processing
#'
#' Data cubes can be processed in parallel where one thread processes one chunk at a time. Setting more threads
#' than the number of chunks of a cube thus has no effect and will not furhter reduce computation times.
#' @param n number of threads
#' @export
gdalcubes_set_threads <- function(n=1) {
  stopifnot(n >= 1)
  stopifnot(n%%1==0)
  libgdalcubes_set_threads(n)
  invisible()
}

#' Query gdalcubes version information
#'
#' @return List with gdalcubes library version information
#' @export
gdalcubes_version <- function() {
  return(libgdalcubes_version())
}

#' Enable or disable debug output from the gdalcubes C++ library
#' @param debug logical, TRUE if you want debug messages
#' @export
gdalcubes_debug_output <- function(debug=TRUE) {
   libgdalcubes_debug_output(debug)
  invisible()
}



#' Enable or disable caching of cubes.
#' 
#' @details
#' Caching has no effect on disk or memory consumption, 
#' it simply tries to reuse existing temporary files where possible.
#' For example, changing only parameters to \code{plot} will not require
#' rerunning the full data cube operation chain.
#' 
#' @param enable logical, TRUE if you want to use the data cube cache
#' @export
gdalcubes_use_cache <- function(enable=TRUE) {
  .pkgenv$use_cube_cache = enable
}

