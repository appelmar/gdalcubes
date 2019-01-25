
#' Create a data cube from an image collection
#' 
#' Create a proxy data cube, which loads data from a given image collection
#'
#' @param image_collection Source image collection as from \code{image_collection} or \code{create_image_collection}
#' @param view A data cube view defining the shape (spatiotemporal extent, resolution, and spatial reference), if missing, a default overview is used
#' @param chunking Vector of length 3 defining the size of data cube chunks in the order time, y, x.
#' @return A proxy data cube object
#' @details 
#' The following steps will be performed when the data cube is requested to read data of a chunk:
#' 
#'  1. Filter images from the input collection that intersect with the spatiotemporal extent of the chunk
#'  2. For all resulting images, apply gdalwarp to reproject / warp the image to the target SRS and size to a GDAL MEM dataset
#'  3. Read the resulting data to the chunk buffer and if pixels already contain non NAN values, apply an aggregation method (as defined in the view) 
#' 
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = cube_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = create_image_collection(L8_files, "L8_L1TP") 
#'  cube(L8.col, v)
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
cube <- function(image_collection, view, chunking=c(16, 256, 256)) {

  stopifnot(is.image_collection(image_collection))
  stopifnot(length(chunking) == 3)
  chunking = as.integer(chunking)
  stopifnot(chunking[1] > 0 && chunking[2] > 0 && chunking[3] > 0)
  
  x = NULL
  if (!missing(view)) {
    stopifnot(is.cube_view(view))
    x = libgdalcubes_create_image_collection_cube(image_collection, as.integer(chunking), view)
  }
  else {
    x = libgdalcubes_create_image_collection_cube(image_collection, as.integer(chunking))
  }
  class(x) <- c("image_collection_cube", "cube", "xptr")
  return(x)
}


is.image_collection_cube <- function(obj) {
  if(!("image_collection_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.cube <- function(obj) {
  if(!("cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}

#' @export
print.cube <- function(x, ...) {
  if (libgdalcubes_is_null(x)) {
    stop("GDAL data cube proxy object is invalid")
  }
  y = libgdalcubes_cube_info(x)
  cat("A GDAL data cube proxy object\n")
  cat("Dimensions:\n")
  print(y$dimensions)
  cat("\n")
  cat("Bands:\n")
  print(y$bands)
  cat("\n")
}

#' Query data cube properties 
#' 
#' @return Size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso dim.cube
#' @param obj a data cube proxy object (class cube)
#' @export
size <- function(obj) {
  if (libgdalcubes_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  x = libgdalcubes_cube_info(obj)
  return(x$size[2:4])
}

#' Query data cube properties 
#' 
#' @return Size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso size
#' @param x a data cube proxy object (class cube)
#' @export
dim.cube <- function(x) {
  return(size(x))
}

#' Query data cube properties 
#' 
#' @return Band names as character vector
#' 
#' @param x a data cube proxy object (class cube)
#' @export
names.cube <- function(x) {
  if (libgdalcubes_is_null(x)) {
    stop("GDAL data cube proxy object is invalid")
  }
  y = libgdalcubes_cube_info(x)
  return(as.character(y$bands$name))
}


#' Query data cube properties 
#' 
#' @return Dimension information as a data.frame, where each row represents a dimension and columns represent properties such as dimension boundaries, names, and chunk size
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
dimensions <- function(obj) {
  if (libgdalcubes_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  x = libgdalcubes_cube_info(obj)
  return(x$dimensions)
}

#' Query data cube properties 
#' 
#' @return A data.frame with rows representing the bands and columns representing properties of a band (name, type, scale, offset, unit)
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
bands <- function(obj) {
  if (libgdalcubes_is_null(obj)) {
    stop("GDAL data cube proxy object is invalid")
  }
  x = libgdalcubes_cube_info(obj)
  return(x$bands)
}

#' Query data cube properties 
#' 
#' @return The spatial reference system expressed as a string readable by GDAL
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
get_projection <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$proj)
}

#' Query data cube properties 
#' 
#' @return Total data size of data cube values expressed in the given unit
#' 
#' @param obj a data cube proxy object (class cube)
#' @param unit Unit of data size, can be "B", "KB", "KiB", "MB", "MiB", "GB", "GiB", "TB", "TiB", "PB", "PiB"
#' @export
get_cubesize <- function(obj, unit="MiB") {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
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
#' @export
nbands <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$size[1])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the time dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
nt <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$size[2])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the y dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
ny <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$size[3])
}

#' Query data cube properties 
#' 
#' @return Number of pixels in the x dimension
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
nx <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$size[4])
}

#' Query data cube properties 
#' 
#' gdalcubes uses a graph (currently a tree) to represent / serialize data cube chains. This function gives a JSON
#' representation, which will be communicated to gdalcubes_server instances to create identical cube instances 
#' remotely.
#' 
#' @return A JSON string representing a graph (currently a tree) that can be used to create the same
#' chain of gdalcubes operations.
#' 
#' @param obj a data cube proxy object (class cube)
#' @export
graph <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$graph)
}
#' 
#' 
#' "cube_view<-" <-function(x,value) {
#'   stopifnot(is.cube(x))
#'   stopifnot(is.cube_view(value))
#'   if (!is.image_collection_cube)
#'     stop("x is no image_collection_cube, updating the data cube view is currently only implemented for image_collection_cube")
#'   libgdalcubes_update_cube_view(x,value)
#'   return(x)
#' }


#' Materialize a data cube as a NetCDF file
#' 
#' This function will read chunks of a data cube and write them to a single NetCDF file.
#' 
#' @param x a data cube proxy object (class cube)
#' @param fname output file name
#' @export
as_ncdf <- function(x, fname = tempfile(pattern = "gdalcubes", fileext = ".nc")) {
  stopifnot(is.cube(x))
  libgdalcubes_eval_cube(x, fname)
}






