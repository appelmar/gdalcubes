
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
#' @param show_progress logical; if TRUE, a progress bar will be shown for actual computations
#' @param default_chunksize length-three vector with chunk size in t, y, x directions or a function taking a data cube size and returning a suggested chunk size 
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
#' gdalcubes_options(threads=1) # reset
#' @export
gdalcubes_options <- function(..., threads, ncdf_compression_level, debug, cache, ncdf_write_bounds, 
                              use_overview_images, show_progress, default_chunksize) {
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
  if (!missing(show_progress)) {
    stopifnot(is.logical(show_progress))
    .pkgenv$show_progress = show_progress
    libgdalcubes_set_progress(show_progress)
  }
  if (!missing(default_chunksize)) {
    if (is.vector(default_chunksize)) {
      stopifnot(length(default_chunksize) == 3)
      stopifnot(all(default_chunksize %% 1 == 0))
    }
    else if (is.function(default_chunksize)) {
      test = default_chunksize(512, 512, 512)
      stopifnot(length(test) == 3)
      stopifnot(all(test %% 1 == 0))
    }
    else {
      stop("Expected a length-three vector or a function for argument default_chunksize")
    }
    .pkgenv$default_chunksize = default_chunksize
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
      use_overview_images = .pkgenv$use_overview_images,
      show_progress = .pkgenv$show_progress,
      default_chunksize = .pkgenv$default_chunksize
    ))
  }
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

#' Set GDAL config options
#' 
#' @details 
#' Details and a list of possible options can be found at 
#' \href{https://gdal.org/user/configoptions.html}{https://gdal.org/user/configoptions.html}.
#' @examples 
#' gdalcubes_set_gdal_config("GDAL_NUM_THREADS", "ALL_CPUS")
#' @export
gdalcubes_set_gdal_config <- function(key, value) {
  stopifnot(length(key) == 1)
  stopifnog(length(value) == 1)
  libgdalcubes_set_gdal_config(as.character(key), as.character(value))
}

#' Calculate a default chunk size based on the cube size and currently used number of threads
#' @examples 
#' .default_chunk_size(12, 1000, 1000)
#' @export
.default_chunk_size <- function(nt, ny, nx) {
  
  nthreads = .pkgenv$threads
  
  ct = 1
  target_pixels_per_chunk = ny * nx * nt / nthreads
  
  # multiples of 256
  cy = max(floor(sqrt(target_pixels_per_chunk) / 256), 1) * 256
  cx = cy
  
  #ar = ny / nx
  #cy = ceiling(sqrt(target_pixels_per_chunk * ar))
  #cx = ceiling(sqrt(target_pixels_per_chunk / ar))
  
  # apply limits
  cx = min(cx, 2048)
  cy = min(cy, 2048)
  cx = max(cx, 256)
  cy = max(cy, 256)
  cx = min(nx, cx)
  cy = min(ny, cy)
  
  return(c(ct, cy, cx))
}



