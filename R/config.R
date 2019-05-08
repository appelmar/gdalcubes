
#' Set the number of threads for parallel data cube processing
#'
#' Data cubes can be processed in parallel where one thread processes one chunk at a time. Setting more threads
#' than the number of chunks of a cube thus has no effect and will not further reduce computation times.
#' @param n number of threads
#' @examples 
#' gdalcubes_set_threads(1)
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
#' @examples 
#' gdalcubes_version()
#' @export
gdalcubes_version <- function() {
  return(libgdalcubes_version())
}

#' Enable or disable debug output from the gdalcubes C++ library
#' @param debug logical, TRUE if you want debug messages
#' @examples 
#' gdalcubes_debug_output(TRUE)
#' gdalcubes_debug_output(FALSE)
#' @export
gdalcubes_debug_output <- function(debug=TRUE) {
   libgdalcubes_debug_output(debug)
  invisible()
}


#' Set compression level for netCDF files produced by gdalcubes
#' @param level integer; compression level, 0 = no compression, 1=fast compression, 9=small compression
#' @examples 
#' gdalcubes_set_ncdf_compression(9) # maximum compression
#' gdalcubes_set_ncdf_compression(0) # no compression
#' @export
gdalcubes_set_ncdf_compression <- function(level=2) {
  stopifnot(level %% 1 == 0)
  stopifnot(level >= 0 && level <= 9)
  .pkgenv$compression_level = level
  invisible()
}

#' Get available GDAL drivers
#' @examples 
#' gdalcubes_gdalformats()
#' @export
gdalcubes_gdalformats <- function() {
  return(libgdalcubes_gdalformats())
}

#' Get the GDAL version used by gdalcubes
#' @examples 
#' gdalcubes_gdalversion()
#' @export
gdalcubes_gdalversion <- function() {
  return(libgdalcubes_gdalversion())
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
#' @examples 
#' gdalcubes_use_cache(FALSE)
#' @export
gdalcubes_use_cache <- function(enable=TRUE) {
  .pkgenv$use_cube_cache = enable
}

