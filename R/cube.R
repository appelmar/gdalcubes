
#' Create a data cube from an image collection
#' 
#' Create a proxy data cube, which loads data from a given image collection according to a data cube view
#'
#' @param image_collection Source image collection as from \code{image_collection} or \code{create_image_collection}
#' @param view A data cube view defining the shape (spatiotemporal extent, resolution, and spatial reference), if missing, a default overview is used
#' @param mask mask pixels of images based on band values, see \code{\link{image_mask}}
#' @param chunking length-3 vector or a function returning a vector of length 3, defining the size of data cube chunks in the order time, y, x.
#' @param incomplete_ok logical, if TRUE (the default), chunks will ignore IO failures and simply use as much images as possible, otherwise the result will contain empty chunks if IO errors or similar occur.
#' @return A proxy data cube object
#' @details 
#' The following steps will be performed when the data cube is requested to read data of a chunk:
#' 
#'  1. Find images from the input collection that intersect with the spatiotemporal extent of the chunk
#'  2. For all resulting images, apply gdalwarp to reproject, resize, and resample to an in-memory GDAL dataset
#'  3. Read the resulting data to the chunk buffer and optionally apply a mask on the result
#'  4. Update pixel-wise aggregator (as defined in the data cube view) to combine values of multiple images within the same data cube pixels
#' 
#' If chunking is provided as a function, it must accept exactly three arguments for the total size of the cube in t, y, and x axes (in this order). 
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' raster_cube(L8.col, v)
#'  
#'  # using a mask on the Landsat quality bit band to filter out clouds
#'  raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
raster_cube <- function(image_collection, view, mask=NULL, chunking=.pkgenv$default_chunksize, incomplete_ok = TRUE) {

  stopifnot(is.image_collection(image_collection))
  if (is.function(chunking)) {
    if (missing(view)) {
      warning("Function to derive chunk sizes is not supported when data cube is missing, using fixed chunk size (1, 512, 512)")
      chunking = c(1, 512, 512)
    }
    else {
      chunking = chunking(view$time$nt, view$space$ny, view$space$nx)
    }
  }
  stopifnot(length(chunking) == 3)
  chunking = as.integer(chunking)
  stopifnot(chunking[1] > 0 && chunking[2] > 0 && chunking[3] > 0)
  if (!is.null(mask)) {
    stopifnot(is.image_mask(mask))
  }
  
  x = NULL
  if (!missing(view)) {
    stopifnot(is.cube_view(view))
    x = gc_create_image_collection_cube(image_collection, as.integer(chunking), mask, !incomplete_ok, view)
  }
  else {
    x = gc_create_image_collection_cube(image_collection, as.integer(chunking), mask, !incomplete_ok)
  }
  class(x) <- c("image_collection_cube", "cube", "xptr")
  return(x)
}



#' Create a data cube from a set of images with the same spatial extent and spatial reference system
#' 
#' Create a spatiotemporal data cube directly from images with identical spatial extent and spatial reference system, similar
#' to a raster stack with an additional dimension supporting both, time and multiple bands / variables.
#'  
#' @details 
#' This function creates a four-dimensional (space, time, bands / variables) raster data cube from a 
#' set of provided files without the need to create an image collection before. This is  possible if all images
#' have the same spatial extent and spatial reference system and can be used for two different file organizations:
#' 
#' 1. If all image files share the same bands / variables, the \code{bands} argument can be ignored (default NULL) can 
#' names of the bands can be specified using the \code{band_names} argument.
#' 
#' 2. If image files represent different band / variable  (e.g. individual files for red, green, and blue channels), the \code{bands}
#' argument must be used to define the corresponding band / variable. Notice that in this case all files are expected to 
#' represent exactly one variable / band at one point in datetime. It is not possible to combine files with different
#' numbers of variables / bands. If image files for different bands have different pixel sizes, the smallest size is used
#' by default.
#' 
#' Notice that to avoid opening all image files in advance,no automatic check whether all images share the 
#' spatial extent and spatial reference system is performed.
#' 
#' @param x character vector where items point to image files
#' @param datetime_values vector of type character, Date, or POSIXct with recording date of images
#' @param bands optional character vector defining the band or spectral band of each item in x, if files relate to different spectral bands or variables
#' @param band_names name of bands, only used if bands is NULL, i.e., if all files contain the same spectral band(s) / variable(s)
#' @param chunking vector of length 3 defining the size of data cube chunks in the order time, y, x.
#' @param dx optional target pixel size in x direction, by default (NULL) the original or highest resolution of images is used 
#' @param dy optional target pixel size in y direction, by default (NULL) the original or highest resolution of images is used 
#' @param incomplete_ok logical, if TRUE (the default), chunks will ignore IO failures and simply use as much images as possible, otherwise the result will contain empty chunks if IO errors or similar occur.
#' @return A proxy data cube object
#' @examples 
#' # toy example, repeating the same image as a daily time series
#' L8_file_nir <- 
#' system.file("L8NY18/LC08_L1TP_014032_20181122_20181129_01_T1/LC08_L1TP_014032_20181122_B5.TIF",
#'             package = "gdalcubes")
#' files = rep(L8_file_nir, 10)
#' datetime = as.Date("2018-11-22") + 1:10
#' stack_cube(files, datetime, band_names = "B05") 
#' 
#' # using a second band from different files
#' L8_file_red <- 
#' system.file("L8NY18/LC08_L1TP_014032_20181122_20181129_01_T1/LC08_L1TP_014032_20181122_B4.TIF",
#'             package = "gdalcubes")
#' files = rep(c(L8_file_nir, L8_file_red), each = 10)
#' datetime = rep(as.Date("2018-11-22") + 1:10, 2)
#' bands = rep(c("B5","B4"), each = 10)
#' stack_cube(files, datetime, bands = bands)         
#'
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.                     
#' @export
stack_cube <- function(x, datetime_values, bands = NULL, band_names = NULL, chunking = c(1, 256, 256), dx=NULL, dy=NULL, incomplete_ok = TRUE) {
  
  
  if (length(datetime_values) != length(x)) {
    stop("x and datetime_values have different length")
  }
  
  if (!is.null(bands)) {
    if (length(bands) != length(x)) {
      stop("x and bands have different length")
    }
  }
  stopifnot(length(chunking) == 3)
  
  if (!is.character(datetime_values)) {
    datetime_values = as.character(datetime_values)
  }
  
  if (!is.null(bands) & !is.null(band_names)) {
    warning("Ignoring band_names because bands have been defined per file")
  }
  
  
  if (is.null(bands)) {
    bands = character(0)
  }
  if (is.null(band_names)) {
    band_names = character(0)
  }
  if (is.null(dx)) {
    dx = -1.0
  }
  if (is.null(dy)) {
    dy = -1.0
  }
 
  x = gc_create_simple_cube(x, datetime_values, bands, band_names, dx, dy, as.integer(chunking), !incomplete_ok)
  class(x) <- c("simple_cube", "cube", "xptr")
  return(x)
}


#' Create a data cube proxy object copy
#' 
#' Copy a data cube proxy object without copying any data
#' 
#' @param  cube source data cube proxy object
#' @return copied data cube proxy object 
#' @details 
#' This internal function copies the complete processing chain / graph of a data cube but does not copy any data
#' It is used internally to avoid in-place modification for operations with potential side effects on source data cubes.
.copy_cube <- function(cube) {
  cc = class(cube)
  cube = gc_copy_cube(cube)
  class(cube) <- cc
  return(cube)
}

#' Read a data cube from a json description file
#' 
#' @param json length-one character vector with a valid json data cube description
#' @param path source data cube proxy object
#' @return data cube proxy object 
#' @details 
#' Data cubes can be stored as JSON description files. These files do not store any data but the recipe
#' how a data cube is constructed, i.e., the chain (or graph) of processes involved. 
#' 
#' Since data cube objects (as returned from \code{\link{raster_cube}}) cannot be saved with normal R methods,
#' the combination of \code{\link{as_json}} and \code{\link{json_cube}} provides a cheap way to save virtual 
#' data cube objects across several R sessions, as in the examples.
#' 
#' @examples
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' cube = raster_cube(L8.col, v) 
#' 
#' # save
#' fname = tempfile()
#' as_json(cube, fname)
#' 
#' # load
#' json_cube(path = fname)  
#' 
#' 
#' @export
json_cube <- function(json, path = NULL) {
  if (!missing(json)) {
    if (!is.null(path)) {
      warning("Expected only one of arguments 'json' and 'path'; path will be ignored")
    }
    cube = gc_from_json_string(json)
  }
  else {
    if (!is.null(path)) {
      cube = gc_from_json_file(path)
    }
    else {
      stop("Missing argument, please provide either a JSON string, or a path to a JSON file")
    }
  }
  class(cube) <- "cube" # TODO: any way to derive exact cube type here?
  return(cube)
}




#' Create a mask for images in a raster data cube 
#'
#' Create an image mask based on a band and provided values to filter pixels of images 
#' read by \code{\link{raster_cube}}
#'
#' @details
#' Values of the selected mask band can be based on a range (by passing \code{min} and \code{max}) or on a set of values (by passing \code{values}). By default
#' pixels with mask values contained in the range or in the values are masked out, i.e. set to NA. Setting \code{invert = TRUE} will invert the masking behavior.
#' Passing \code{values} will override \code{min} and \code{max}.
#' 
#' @note 
#' Notice that masks are applied per image while reading images as a raster cube. They can be useful to eliminate e.g. cloudy pixels before applying the temporal aggregation to
#' merge multiple values for the same data cube pixel.
#' 
#' @examples 
#' image_mask("SCL", values = c(3,8,9)) # Sentinel 2 L2A: mask cloud and cloud shadows
#' image_mask("BQA", bits=4, values=16) # Landsat 8: mask clouds
#' image_mask("B10", min = 8000, max=65000) 
#' 
#' @param band name of the mask band
#' @param min minimum value, values between \code{min} and \code{max} will be masked
#' @param max maximum value, values between \code{min} and \code{max} will be masked 
#' @param values numeric vector; specific values that will be masked. 
#' @param bits for bitmasks, extract the given bits (integer vector) with a bitwise AND before filtering the mask values, bit indexes are zero-based
#' @param invert logical; invert mask
#' @export
image_mask <- function(band, min=NULL, max=NULL, values=NULL, bits=NULL, invert=FALSE) {
  if (is.null(values) && is.null(min) && is.null(max)) {
    stop("either values or min and max must be provided")
  } 
  if (is.null(values) && is.null(min) && !is.null(max)) {
    stop("either values or min AND max must be provided")
  } 
  if (is.null(values) && !is.null(min) && is.null(max)) {
    stop("either values or min AND max must be provided")
  } 
  if (!is.null(values)) {
    if (!is.null(min) || !is.null(max)) {
      warning("using values instead of min / max")
    }
    out = list(band = band, values = values, invert = invert, bits = bits)
  }  
  else {
    out = list(band = band, min = min, max = max, invert = invert, bits = bits)
  }
  class(out) <- "image_mask"
  return(out)
}

is.image_mask <- function(obj) {
  if(!("image_mask" %in% class(obj))) {
    return(FALSE)
  }
  return(TRUE)
}


is.image_collection_cube <- function(obj) {
  if(!("image_collection_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.cube <- function(obj) {
  if(!("cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}

#' Print data cube information
#' 
#' Prints information about the dimensions and bands of a data cube.
#' 
#' @param x Object of class "cube"
#' @param ... Further arguments passed to the generic print function
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' print(raster_cube(L8.col, v))
#' @export
print.cube <- function(x, ...) {
  if (gc_is_null(x)) {
    stop("Data cube proxy object is invalid")
  }
  y = gc_cube_info(x)
  cat("A data cube proxy object\n")
  cat("\n")
  cat("Dimensions:\n")
  dimensions = data.frame(
    #name = c("time","y","x"),
    low = sapply(y$dimensions, function(z) z$low),
    high = sapply(y$dimensions, function(z) z$high),
    count = sapply(y$dimensions, function(z) z$count),
    pixel_size = sapply(y$dimensions, function(z) z$pixel_size),
    chunk_size = sapply(y$dimensions, function(z) z$chunk_size)
  )
  if (!is.null(y$dimensions$t$values)) {
    nmax = 5
    str = paste(head(y$dimensions$t$values,nmax), collapse=",")
    if (length(y$dimensions$t$values) > nmax)
      str = paste0(str, ",...")
    dimensions$values = c(str, "","")
  }
  rownames(dimensions) = c("t","y","x")
  print(dimensions)
  
  cat("\n")
  cat("Bands:\n")
  print(y$bands)
  cat("\n")
}

#' Query data cube properties 
#' 
#' @return size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso \code{\link{dim.cube}}
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' size(raster_cube(L8.col, v))
#' @export
size <- function(obj) {
  if (gc_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  x = gc_cube_info(obj)
  return(x$size[2:4])
}

#' Query data cube properties 
#' 
#' @return size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso \code{\link{size}}
#' @param x a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' dim(raster_cube(L8.col, v))
#' @export
dim.cube <- function(x) {
  return(size(x))
}

#' Query data cube properties 
#' 
#' @return Band names as character vector
#' 
#' @param x a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' names(raster_cube(L8.col, v))
#' @export
names.cube <- function(x) {
  if (gc_is_null(x)) {
    stop("GDAL data cube proxy object is invalid")
  }
  y = gc_cube_info(x)
  return(as.character(y$bands$name))
}


#' Query data cube properties 
#' 
#' @return Dimension information as a list
#' 
#' @details Elements of the returned list represent individual dimensions with properties such as dimension boundaries, names, and chunk size stored as inner lists
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' dimensions(raster_cube(L8.col, v))
#' @export
dimensions <- function(obj) {
  if (gc_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  y = gc_cube_info(obj)
  return(y$dimensions)
}

#' Query data cube properties 
#' 
#' @return A data.frame with rows representing the bands and columns representing properties of a band (name, type, scale, offset, unit)
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' bands(raster_cube(L8.col, v))
#' @export
bands <- function(obj) {
  if (gc_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  x = gc_cube_info(obj)
  return(x$bands)
}

#' Query data cube properties 
#' 
#' @return The spatial reference system expressed as a string readable by GDAL
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' srs(raster_cube(L8.col, v))
#' @export
srs <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$srs)
}

#' Query data cube properties 
#' 
#' @return The spatial reference system expressed as proj4 string
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' proj4(raster_cube(L8.col, v))
#' @export
proj4 <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$proj4)
}


#' Query data cube properties 
#' 
#' @return Total data size of data cube values expressed in the given unit
#' 
#' @param obj a data cube proxy object (class cube)
#' @param unit Unit of data size, can be "B", "KB", "KiB", "MB", "MiB", "GB", "GiB", "TB", "TiB", "PB", "PiB"
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' memsize(raster_cube(L8.col, v))
#' @export
memsize <- function(obj, unit="MiB") {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  size_bytes = prod(x$size) * 8 # assuming everything is double
  return(switch(unit,
         B = size_bytes,
         KB = size_bytes / 1000,
         KiB = size_bytes / 1024,
         MB = size_bytes / (1000^2),
         MiB = size_bytes / (1024^2),
         GB = size_bytes / (1000^3),
         GiB = size_bytes / (1024^3),
         TB = size_bytes / (1000^4),
         TiB = size_bytes / (1024^4),
         PB = size_bytes / (1000^5),
         PiB = size_bytes / (1024^5)
  ))
}


#' Query data cube properties 
#' 
#' @return Number of bands
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' nbands(raster_cube(L8.col, v))
#' @export
nbands <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$size[1])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the time dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' nt(raster_cube(L8.col, v))
#' @export
nt <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$size[2])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the y dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' ny(raster_cube(L8.col, v))
#' @export
ny <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$size[3])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the x dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' nx(raster_cube(L8.col, v))
#' @export
nx <- function(obj) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  return(x$size[4])
}

#' Query data cube properties 
#' 
#' gdalcubes internally uses a graph to serialize data cubes (including chained operations on cubes). This function derives a JSON
#' representation, which can be used to save data cube objects without pixel data to disk.
#' 
#' @return If file is NULL, the function returns a JSON string representing a graph that can be used to recreate the same
#' chain of gdalcubes operations even in a different R sessions.
#' 
#' @param obj a data cube proxy object (class cube)
#' @param file optional output file 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' cat(as_json(select_bands(raster_cube(L8.col, v), c("B04", "B05"))))
#' @export
as_json <- function(obj, file = NULL) {
  stopifnot(is.cube(obj))
  x = gc_cube_info(obj)
  g = jsonlite::prettify(x$graph)
  if (!is.null(file)) {
    writeLines(g, file)
    return(invisible())
  }
  return(g)
}

      


#' Helper function to define packed data exports by min / max values 
#'
#' This function can be used to define packed exports in \code{\link{write_ncdf}}
#' and \code{\link{write_tif}}. It will generate scale and offset values with maximum precision (unless simplify=TRUE).
#'
#' @details
#' Nodata values will be mapped to the lowest value of the target data type.
#' 
#' Arguments min and max must have length 1 or length equal to the number of bands of the data cube to be exported. In the former
#' case, the same values are used for all bands of the exported target cube, whereas the latter case allows to use different 
#' ranges for different bands.
#' 
#' @note 
#' Using simplify=TRUE will round scale values to the next smaller power of 10.
#' 
#' @examples 
#' ndvi_packing = pack_minmax(type="int16", min=-1, max=1)
#' ndvi_packing
#' 
#' @param type target data type of packed values (one of "uint8", "uint16", "uint32", "int16", or "int32")
#' @param min numeric; minimum value(s) of original values, will be packed to the 2nd lowest value of the target data type
#' @param max numeric; maximum value(s) in original scale, will be packed to the highest value of the target data type
#' @param simplify logical; round resulting scale and offset to power of 10 values
#' @export
pack_minmax <- function(type="int16", min, max, simplify=FALSE) {
  
  stopifnot(length(min) == length(max))
  
  if (type == "int16") {
    nodata = -2^15
    low = -2^15+1
    high = 2^15 - 1
    scale = (max-min)/(high-low)
    offset = min - low * scale 
    out = list(type="int16", offset=offset, scale=scale, nodata=nodata)
  }
  else if (type == "int32") {
    nodata = -2^31
    low = -2^(31)+1
    high = 2^31 - 1
    scale = (max-min)/(high-low)
    offset = min - low * scale 
    out = list(type="int32", offset=offset, scale=scale, nodata=nodata)
  }
  else if (type == "uint8") {
    nodata = 0
    low = 1
    high = 2^8 - 1
    scale = (max-min)/(high-low)
    offset = min - low * scale 
    out = list(type="uint8", offset=offset, scale=scale, nodata=nodata)
  }
  else if (type == "uint16") {
    nodata = 0
    low =  1
    high = 2^16 - 1
    scale = (max-min)/(high-low)
    offset = min - low * scale 
    out = list(type="uint16", offset=offset, scale=scale, nodata=nodata)
  }
  else if (type == "uint32") {
    nodata =  0
    low =  1
    high = 2^32 - 1
    scale = (max-min)/(high-low)
    offset = min - low * scale 
    out = list(type="uint32", offset=offset, scale=scale, nodata=nodata)
  }
  else {
    stop("Invalid data type for packed export.")
  }
  
  if (simplify) {
    floor_10 <- function(x) 10^floor(log10(x))
    out$scale = floor_10(out$scale)
  }
  return(out)
}




#' Export a data cube as netCDF file(s)
#' 
#' This function will read chunks of a data cube and write them to a single (the default) or multitple (if \code{chunked = TRUE}) netCDF file(s). The resulting
#' file(s) uses the enhanced netCDF-4 format, supporting chunking and compression.
#' 
#' @seealso \code{\link{gdalcubes_options}} 
#' @param x a data cube proxy object (class cube)
#' @param fname output file name
#' @param overwrite logical; overwrite output file if it already exists
#' @param write_json_descr logical; write a JSON description of x as additional file
#' @param with_VRT logical; write additional VRT datasets (one per time slice)
#' @param pack reduce output file size by packing values (see Details), defaults to no packing
#' @param chunked logical; if TRUE, write one netCDF file per chunk; defaults to FALSE 
#' 
#' @seealso \code{\link{pack_minmax}}
#' 
#' @details 
#' The resulting netCDF file(s) contain three dimensions (t, y, x) and bands as variables.
#' 
#' If \code{write_json_descr} is TRUE, the function will write an addition file with the same name as the NetCDF file but 
#' ".json" suffix. This file includes a serialized description of the input data cube, including all chained data cube operations.
#'
#' To reduce the size of created files, values can be packed by applying a scale factor and an offset value and using a smaller
#' integer data type for storage (only supported if \code{chunked = TRUE}). The \code{pack} argument can be either NULL (the default), or a list with elements \code{type}, \code{scale}, \code{offset}, 
#' and \code{nodata}. \code{type} can be any of "uint8", "uint16" , "uint32", "int16", or "int32". \code{scale}, \code{offset}, and 
#' \code{nodata} must be numeric vectors with length one or length equal to the number of data cube bands (to use different values for different bands). 
#' The helper function  \code{\link{pack_minmax}} can be used to derive offset and scale values with maximum precision from minimum and maximum data values on
#' original scale.
#' 
#' If \code{chunked = TRUE}, names of the produced files will start with \code{name} (with removed extension), followed by an underscore and the internal integer chunk number. 
#' 
#' @note Packing is currently ignored if \code{chunked = TRUE}
#' 
#' @return returns (invisibly) the path of the created netCDF file(s) 
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' write_ncdf(select_bands(raster_cube(L8.col, v), c("B04", "B05")), fname=tempfile(fileext = ".nc"))
#' @export
write_ncdf <- function(x, fname = tempfile(pattern = "gdalcubes", fileext = ".nc"), overwrite = FALSE, 
                       write_json_descr = FALSE, with_VRT = FALSE, pack = NULL, chunked = FALSE) {
  stopifnot(is.cube(x))
  fname = path.expand(fname)
  if (!overwrite && file.exists(fname)) {
    stop("File already exists, please change the output filename or set overwrite = TRUE")
  }
  
  
  if (!is.null(pack)) {
    stopifnot(is.list(pack))
    stopifnot(length(pack$offset) == 1 || length(pack$offset) == nbands(x))
    stopifnot(length(pack$scale) == 1 || length(pack$scale) == nbands(x))
    stopifnot(length(pack$nodata) == 1 || length(pack$nodata) == nbands(x))
    stopifnot(length(pack$offset) == length(pack$scale))
    stopifnot(length(pack$offset) == length(pack$nodata))
  }
  
  if (!is.null(pack) && chunked) {
    warning("Since chunked = TRUE, packing will be ignored (data type will remain 8 byte double)")
  }
  if (.pkgenv$ncdf_write_bounds && chunked) {
    warning("Since chunked = TRUE, resulting netCDF files will not include bounds variables.")
  }
  
  if (!chunked) {
    if (.pkgenv$use_cube_cache) {
      j = gc_simple_hash(as_json(x))
      if (!is.null(.pkgenv$cube_cache[[j]])
          && file.exists(.pkgenv$cube_cache[[j]])) {
        file.copy(from=.pkgenv$cube_cache[[j]], to = fname, overwrite=TRUE)
      }
      else {
        gc_eval_cube(x, fname, .pkgenv$compression_level, with_VRT, .pkgenv$ncdf_write_bounds, pack)
      }
    }
    else {
      gc_eval_cube(x, fname, .pkgenv$compression_level, with_VRT, .pkgenv$ncdf_write_bounds, pack)
    }
  }
  else {
    gc_write_chunks_ncdf(x, dirname(fname), tools::file_path_sans_ext(basename(fname)), .pkgenv$compression_level)
  }
  
  if (write_json_descr) {
    writeLines(as_json(x), paste(fname, ".json", sep=""))
  }
  if (!chunked) {
    invisible(fname)
  }
  else {
    list.files(dirname(fname), pattern=paste(tools::file_path_sans_ext(basename(fname)), "_[0-9]+.nc", sep=""), full.names = TRUE)
  }
}







#' Export a data cube as a collection of GeoTIFF files
#' 
#' This function will time slices of a data cube as GeoTIFF files
#' in a given directory. 
#' 
#' @param x a data cube proxy object (class cube)
#' @param dir destination directory
#' @param prefix output file name
#' @param overviews logical; generate overview images 
#' @param COG logical; create cloud-optimized GeoTIFF files (forces overviews=TRUE)
#' @param rsmpl_overview resampling method for overviews (image pyramid) generation (see \url{https://gdal.org/en/stable/programs/gdaladdo.html} for available methods)
#' @param creation_options additional creation options for resulting GeoTIFF files, e.g. to define compression (see \url{https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options})
#' @param write_json_descr logical; write a JSON description of x as additional file
#' @param pack reduce output file size by packing values (see Details), defaults to no packing
#' 
#' @seealso \code{\link{pack_minmax}}
#' 
#' @return  returns (invisibly) a vector of paths pointing to the created GeoTIFF files
#' 
#' @details 
#' 
#' If \code{write_json_descr} is TRUE, the function will write an additional file with name according to prefix (if not missing) or simply cube.json 
#' This file includes a serialized description of the input data cube, including all chained data cube operations.
#' 
#' Additional GDAL creation options for resulting GeoTIFF files must be passed as a named list of simple strings, where element names refer to the key. For example,
#' \code{creation_options = list("COMPRESS" = "DEFLATE", "ZLEVEL" = "5")} would enable deflate compression at level 5.
#' 
#' To reduce the size of created files, values can be packed by applying a scale factor and an offset value and using a smaller
#' integer data type for storage. The \code{pack} argument can be either NULL (the default), or a list with elements \code{type}, \code{scale}, \code{offset}, 
#' and \code{nodata}. \code{type} can be any of "uint8", "uint16" , "uint32", "int16", or "int32". \code{scale}, \code{offset}, and 
#' \code{nodata} must be numeric vectors with length one or length equal to the number of data cube bands (to use different values for different bands). 
#' The helper function  \code{\link{pack_minmax}} can be used to derive offset and scale values with maximum precision from minimum and maximum data values on
#' original scale.
#' 
#' If \code{overviews=TRUE}, the numbers of pixels are halved until the longer spatial dimensions counts less than 256 pixels.
#' Setting \code{COG=TRUE} automatically sets \code{overviews=TRUE}.
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' write_tif(select_bands(raster_cube(L8.col, v), c("B04", "B05")), dir=tempdir())
#' @export
write_tif <- function(x, dir = tempfile(pattern=""), prefix = basename(tempfile(pattern = "cube_")), overviews = FALSE, 
                      COG = FALSE, rsmpl_overview="nearest", creation_options = NULL , write_json_descr=FALSE, pack = NULL) {
  stopifnot(is.cube(x))
  dir = path.expand(dir)
  # if (dir.exists(dir) && prefix == "" && length(list.files(dir, include.dirs = TRUE) > 0)) {
  #   stop("Directory already exists and is not empty, please either")
  # }
  
  if (!(is.null(creation_options) || is.list(creation_options))) {
    stop("Expected either NULL or a list as creation_options argument.")
  }
  
  if (!is.character(rsmpl_overview)) {
    stop("Expected a chracte as rsmpl_overview argument.")
  }
  
  if (!overviews && COG) {
    overviews = TRUE
  }
  
  if (!is.null(pack)) {
    stopifnot(is.list(pack))
    stopifnot(length(pack$offset) == 1 || length(pack$offset) == nbands(x))
    stopifnot(length(pack$scale) == 1 || length(pack$scale) == nbands(x))
    stopifnot(length(pack$nodata) == 1 || length(pack$nodata) == nbands(x))
    stopifnot(length(pack$offset) == length(pack$scale))
    stopifnot(length(pack$offset) == length(pack$nodata))
  }
  
  
  # TODO: find out how to enable caching
  gc_write_tif(x, dir, prefix, overviews, COG, creation_options, rsmpl_overview,  pack)
  if (write_json_descr) {
    if (prefix == "") {
      writeLines(as_json(x), file.path(dir, "cube.json"))
    }
    else {
      writeLines(as_json(x), file.path(dir, paste(prefix, ".json", sep="")))
    }
  }
  return(invisible(list.files(path = dir,pattern = paste(prefix, ".*\\.tif", sep=""), full.names = TRUE)))
}













#' Query coordinate values for all dimensions of a data cube 
#' 
#' Dimension values give the coordinates along the spatial and temporal axes of a data cube.
#' 
#' @param obj a data cube proxy (class cube), or a data cube view object
#' @param datetime_unit unit used to format values in the datetime dimension, one of "Y", "m", "d", "H", "M", "S", defaults to the unit of the cube.
#' @return list with elements t,y,x
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' dimension_values(raster_cube(L8.col, v))
#' @export
dimension_values <- function(obj, datetime_unit=NULL) {
  if (is.cube(obj)) {
    if (is.null(datetime_unit)) {
      datetime_unit = ""
    }
    return(gc_dimension_values(obj, datetime_unit)) 
  }
  else if (is.cube_view(obj)) {
    if (is.null(datetime_unit)) {
      datetime_unit = ""
    }    
    return(gc_dimension_values_from_view(obj, datetime_unit)) 
  }
  else {
    stop("obj must be either from class cube or from class cube_view")
  }
}



#' Query coordinate bounds for all dimensions of a data cube 
#' 
#' Dimension values give the coordinates bounds the spatial and temporal axes of a data cube.
#' 
#' @param obj a data cube proxy (class cube)
#' @param datetime_unit unit used to format values in the datetime dimension, one of "Y", "m", "d", "H", "M", "S", defaults to the unit of the cube.
#' @return list with elements t,y,x, each a list with two elements, start and end
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' dimension_bounds(raster_cube(L8.col, v))
#' @export
dimension_bounds <- function(obj, datetime_unit=NULL) {
  stopifnot(is.cube(obj))
  if (is.null(datetime_unit)) {
    datetime_unit = ""
  }
  bnds = gc_dimension_bounds(obj, datetime_unit)
  out = list(t = list(start = bnds$t[seq(1,length(bnds$t), by = 2)], end = bnds$t[seq(2,length(bnds$t), by = 2)]),
             y = list(start = bnds$y[seq(1,length(bnds$y), by = 2)], end = bnds$y[seq(2,length(bnds$y), by = 2)]),
             x = list(start = bnds$x[seq(1,length(bnds$x), by = 2)], end = bnds$x[seq(2,length(bnds$x), by = 2)]))
  return(out) 
}

