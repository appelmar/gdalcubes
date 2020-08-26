
#' Set or read global options of the gdalcubes package
#'
#' Set global package options to change the default behavior of gdalcubes. These include how many threads are used to process data cubes, how created netCDF files are compressed, and whether
#' or not debug messages should be printed.
#'
#' @param ... not used
#' @param threads number of threads used to process data cubes
#' @param ncdf_compression_level integer; compression level for created netCDF files, 0=no compression, 1=fast compression, 9=small compression
#' @param debug logical;  print debug messages
#' @param cache logical; TRUE if temporary data cubes should be cached to support fast reprocessing of the same cubes
#' @param ncdf_write_bounds logical; write dimension bounds as additional variables in netCDF files
#' @param use_overview_images logical; if FALSE, all images are read on original resolution and existing overviews will be ignored
#' @details 
#' Data cubes can be processed in parallel where one thread processes one chunk at a time. Setting more threads
#' than the number of chunks of a cube thus has no effect and will not further reduce computation times.
#' 
#' Caching has no effect on disk or memory consumption, 
#' it simply tries to reuse existing temporary files where possible.
#' For example, changing only parameters to \code{plot} will not require
#' rerunning the full data cube operation chain.
#' 
#' Passing no arguments will return the current options as a list.
#' @examples 
#' gdalcubes_options(threads=4) # set the number of threads
#' gdalcubes_options() # print current options
#' @export
gdalcubes_options <- function(..., threads, ncdf_compression_level, debug, cache, ncdf_write_bounds, use_overview_images) {
  if (!missing(threads)) {
    stopifnot(threads >= 1)
    stopifnot(threads%%1==0)
    libgdalcubes_set_threads(threads)
    .pkgenv$threads = threads
  }
  if (!missing(ncdf_compression_level)) {
    stopifnot(ncdf_compression_level %% 1 == 0)
    stopifnot(ncdf_compression_level >= 0 && ncdf_compression_level <= 9)
    .pkgenv$compression_level = ncdf_compression_level
  }
  if (!missing(debug)) {
    stopifnot(is.logical(debug))
    libgdalcubes_debug_output(debug)
    .pkgenv$debug = debug
  }
  if (!missing(cache)) {
    stopifnot(is.logical(cache))
    .pkgenv$use_cube_cache = cache
  }
  if (!missing(ncdf_write_bounds)) {
    stopifnot(is.logical(ncdf_write_bounds))
    .pkgenv$ncdf_write_bounds = ncdf_write_bounds
  }
  if (!missing(use_overview_images)) {
    stopifnot(is.logical(use_overview_images))
    .pkgenv$use_overview_images = use_overview_images
    libgdalcubes_set_use_overviews(use_overview_images)
  }

  
  # if (!missing(swarm)) {
  #   stopifnot(is.character(swarm))
  #   # check whether all endpoints are accessible
  #   #libgdalcubes_set_swarm(swarm)
  #   warning("swarm mode is currently not supported by the R package")
  #   .pkgenv$swarm = swarm
  # }
  if (nargs() == 0) {
    return(list(
      threads = .pkgenv$threads,
      ncdf_compression_level = .pkgenv$compression_level,
      debug = .pkgenv$debug,
      cache = .pkgenv$use_cube_cache,
      ncdf_write_bounds = .pkgenv$ncdf_write_bounds,
      use_overview_images = .pkgenv$use_overview_images
    ))
  }
}



#' Set the number of threads for parallel data cube processing
#'
#' Data cubes can be processed in parallel where one thread processes one chunk at a time. Setting more threads
#' than the number of chunks of a cube thus has no effect and will not further reduce computation times.
#' @param n number of threads
#' @examples 
#' gdalcubes_set_threads(1)
#' @note THIS FUNCTION IS DEPRECATED AND IS GOING TO BE REPLACED BY \code{\link{gdalcubes_options}}.
#' @export
gdalcubes_set_threads <- function(n=1) {
  stopifnot(n >= 1)
  stopifnot(n%%1==0)
  libgdalcubes_set_threads(n)
  .pkgenv$threads = n
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
#' @note THIS FUNCTION IS DEPRECATED AND IS GOING TO BE REPLACED BY \code{\link{gdalcubes_options}}.
#' @export
gdalcubes_debug_output <- function(debug=TRUE) {
   libgdalcubes_debug_output(debug)
  .pkgenv$debug = debug
  invisible()
}


#' Set compression level for netCDF files produced by gdalcubes
#' @param level integer; compression level, 0 = no compression, 1=fast compression, 9=small compression
#' @examples 
#' gdalcubes_set_ncdf_compression(9) # maximum compression
#' gdalcubes_set_ncdf_compression(0) # no compression
#' @note THIS FUNCTION IS DEPRECATED AND IS GOING TO BE REPLACED BY \code{\link{gdalcubes_options}}.
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

#' Check if GDAL was built with GEOS
#' @examples 
#' gdalcubes_gdal_has_geos()
#' @export
gdalcubes_gdal_has_geos <- function() {
  return(libgdalcubes_gdal_has_geos())
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
#' @note THIS FUNCTION IS DEPRECATED AND IS GOING TO BE REPLACED BY \code{\link{gdalcubes_options}}.
#' @export
gdalcubes_use_cache <- function(enable=TRUE) {
  .pkgenv$use_cube_cache = enable
}

