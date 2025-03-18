
#' Set or read global options of the gdalcubes package
#'
#' Set global package options to change the default behavior of gdalcubes. These include how many parallel processes are used
#' to process data cubes, how created netCDF files are compressed, and whether or not debug messages should be printed.
#'
#' @param ... not used
#' @param parallel number of parallel workers used to process data cubes or TRUE to use the number of available cores automatically
#' @param ncdf_compression_level integer; compression level for created netCDF files, 0=no compression, 1=fast compression, 9=small compression
#' @param debug logical;  print debug messages
#' @param cache logical; TRUE if temporary data cubes should be cached to support fast reprocessing of the same cubes
#' @param ncdf_write_bounds logical; write dimension bounds as additional variables in netCDF files
#' @param use_overview_images logical; if FALSE, all images are read on original resolution and existing overviews will be ignored
#' @param show_progress logical; if TRUE, a progress bar will be shown for actual computations
#' @param default_chunksize length-three vector with chunk size in t, y, x directions or a function taking a data cube size and returning a suggested chunk size 
#' @param streaming_dir directory where temporary binary files for process streaming will be written to
#' @param log_file character, if empty string or NULL, diagnostic messages will be printed to the console, otherwise to the provided file
#' @param threads number of threads used to process data cubes (deprecated)
#' @details 
#' Data cubes can be processed in parallel where the number of chunks in a cube is distributed among parallel
#' worker processes. The actual number of used workers can be lower if a data cube as less chunks. If parallel
#' is TRUE, the number of available cores is used. Setting parallel = FALSE can be used to disable parallel processing.
#' Notice that since version 0.6.0, separate processes are being used instead of parallel threads to avoid 
#' possible R session crashes due to some multithreading issues. 
#' 
#' Caching has no effect on disk or memory consumption, 
#' it simply tries to reuse existing temporary files where possible.
#' For example, changing only parameters to \code{plot} will void
#' reprocessing the same data cube if cache is TRUE.
#' 
#' The streaming directory can be used to control the performance of user-defined functions,
#' if disk IO is a bottleneck. Ideally, this can be set to a directory on a shared memory device.
#' 
#' Passing no arguments will return the current options as a list.
#' @examples 
#' gdalcubes_options(parallel=4) # set the number 
#' gdalcubes_options() # print current options
#' gdalcubes_options(parallel=FALSE) # reset
#' @export
gdalcubes_options <- function(..., parallel, ncdf_compression_level, debug, cache, ncdf_write_bounds, 
                              use_overview_images, show_progress, default_chunksize, streaming_dir, 
                              log_file, threads) {
  if (!missing(threads)) {
    .Deprecated("parallel","gdalcubes", "'threads' option is deprecated; please use 'parallel' instead")
    parallel = threads
  }
  if (!missing(parallel)) {
    if (is.logical(parallel)) {
      if (!parallel) {
        parallel = 1
      }
      else {
        parallel = gc_detect_cores()
        if (parallel == 0) {
          warning("Could not detect the number of available cores automaticall, please set manually")
          parallel = .pkgenv$parallel # use current value
        }
      }
    }
    stopifnot(parallel >= 1)
    stopifnot(parallel%%1==0)
    .pkgenv$parallel = parallel
   
    gc_set_process_execution(.pkgenv$parallel, .pkgenv$worker.cmd, .pkgenv$worker.debug, .pkgenv$worker.compression_level, 
                             .pkgenv$worker.use_overview_images, .pkgenv$worker.gdal_options)
  
  }
  if (!missing(ncdf_compression_level)) {
    stopifnot(ncdf_compression_level %% 1 == 0)
    stopifnot(ncdf_compression_level >= 0 && ncdf_compression_level <= 9)
    .pkgenv$compression_level = ncdf_compression_level
    # This option does NOT automatically set worker.compression_level
  }
  if (!missing(debug)) {
    stopifnot(is.logical(debug))
    .pkgenv$debug = debug
    .pkgenv$worker.debug = debug # set debug mode for worker processes, too
   
    gc_set_process_execution(.pkgenv$parallel, .pkgenv$worker.cmd, .pkgenv$worker.debug, .pkgenv$worker.compression_level, 
                             .pkgenv$worker.use_overview_images, .pkgenv$worker.gdal_options)
    
    gc_set_err_handler(.pkgenv$debug, .pkgenv$log_file)
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
    .pkgenv$worker.use_overview_images = use_overview_images # set value for worker processes, too
   
    gc_set_process_execution(.pkgenv$parallel, .pkgenv$worker.cmd, .pkgenv$worker.debug, .pkgenv$worker.compression_level, 
                             .pkgenv$worker.use_overview_images, .pkgenv$worker.gdal_options)
  
    gc_set_use_overviews(use_overview_images)
  }
  if (!missing(show_progress)) {
    stopifnot(is.logical(show_progress))
    .pkgenv$show_progress = show_progress
    gc_set_progress(show_progress)
  }
  if (!missing(streaming_dir)) {
    stopifnot(is.character(streaming_dir))
    .pkgenv$streaming_dir = streaming_dir
    gc_set_streamining_dir(streaming_dir)
  }
  if (!missing(log_file)) {
    if (is.null(log_file)) log_file = ""
    if (is.na(log_file)) log_file = ""
    stopifnot(is.character(log_file))
    stopifnot(length(log_file) == 1)
    
    if (nchar(log_file[1]) > 0) {
      if (dir.exists(log_file[1])) {
        stop("provided log file is a directory")
      }
      if (!dir.exists(dirname(log_file[1]))) {
        stop("parent directory of provided log file does not exist")
      }
      if (endsWith(log_file[1],.Platform$file.sep)) {
        stop("provided log file is a directory")
      }
    }
    .pkgenv$log_file = log_file
    gc_set_err_handler(.pkgenv$debug, .pkgenv$log_file)
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
  #   #gc_set_swarm(swarm)
  #   warning("swarm mode is currently not supported by the R package")
  #   .pkgenv$swarm = swarm
  # }
  if (nargs() == 0) {
    return(list(
      parallel = .pkgenv$parallel,
      ncdf_compression_level = .pkgenv$compression_level,
      debug = .pkgenv$debug,
      cache = .pkgenv$use_cube_cache,
      ncdf_write_bounds = .pkgenv$ncdf_write_bounds,
      use_overview_images = .pkgenv$use_overview_images,
      show_progress = .pkgenv$show_progress,
      default_chunksize = .pkgenv$default_chunksize,
      streaming_dir = .pkgenv$streaming_dir
    ))
  }
}

#' Get available GDAL drivers
#' @examples 
#' gdalcubes_gdalformats()
#' @export 
gdalcubes_gdalformats <- function() {
  return(gc_gdalformats())
}

#' Check if GDAL was built with GEOS
#' @examples 
#' gdalcubes_gdal_has_geos()
#' @export
gdalcubes_gdal_has_geos <- function() {
  return(gc_gdal_has_geos())
}


#' Get the GDAL version used by gdalcubes
#' @examples 
#' gdalcubes_gdalversion()
#' @export
gdalcubes_gdalversion <- function() {
  return(gc_gdalversion())
}

#' Set GDAL config options
#' 
#' @param key name of a GDAL config option to be set
#' @param value value
#' @details 
#' Details and a list of possible options can be found at 
#' \href{https://gdal.org/en/stable/user/configoptions.html}{https://gdal.org/en/stable/user/configoptions.html}.
#' @examples 
#' gdalcubes_set_gdal_config("GDAL_NUM_THREADS", "ALL_CPUS")
#' @export
gdalcubes_set_gdal_config <- function(key, value) {
  # TODO: implement ... as named elements to set several GDAL options with one function call
  stopifnot(length(key) == 1)
  stopifnot(length(value) == 1)
  .pkgenv$worker.gdal_options[as.character(key)] = as.character(value)
  gc_set_gdal_config(as.character(key), as.character(value))
  gc_set_process_execution(.pkgenv$parallel, .pkgenv$worker.cmd, .pkgenv$worker.debug, .pkgenv$worker.compression_level, 
                           .pkgenv$worker.use_overview_images, .pkgenv$worker.gdal_options)

}

#' Calculate a default chunk size based on the cube size and currently used number of thread
#' @param nt size of a cube in time direction
#' @param ny size of a cube in y direction
#' @param nx size of a cube in x direction
#' @examples 
#' .default_chunk_size(12, 1000, 1000)
#' @noRd
.default_chunk_size <- function(nt, ny, nx) {
  nparallel = .pkgenv$parallel
  target_nchunks_space = ceiling(2*nparallel)
  cx = sqrt(nx*ny / target_nchunks_space)
  cx = ceiling(cx / 64)*64 # use multiples of 64 only
  cy = cx
  
  # apply limits
  cx = min(cx, 1024)
  cy = min(cy, 1024)
  cx = max(cx, 64)
  cy = max(cy, 64)
  cx = min(nx, cx)
  cy = min(ny, cy)
  
  return(c(1, ceiling(cy), ceiling(cx)))
}



