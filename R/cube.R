
#' Create a data cube from an image collection
#' 
#' Create a proxy data cube, which loads data from a given image collection according to a data cube view
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
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' raster_cube(L8.col, v)
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
raster_cube <- function(image_collection, view, chunking=c(16, 256, 256)) {

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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' print(raster_cube(L8.col, v))
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
#' @return size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso \code{\link{dim.cube}}
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' size(raster_cube(L8.col, v))
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
#' @return size of a data cube (number of cells) as integer vector in the order t, y, x
#' @seealso \code{\link{size}}
#' @param x a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' names(raster_cube(L8.col, v))
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
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' dimensions(raster_cube(L8.col, v))
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
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' bands(raster_cube(L8.col, v))
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
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
  x = libgdalcubes_cube_info(obj)
  return(x$size[4])
}

#' Query data cube properties 
#' 
#' gdalcubes uses a graph (currently a tree) to serialize data cubes (including chains of cubes). This function gives a JSON
#' representation, which will be communicated to gdalcubes_server instances to create identical cube instances 
#' remotely.
#' 
#' @return A JSON string representing a graph (currently a tree) that can be used to create the same
#' chain of gdalcubes operations.
#' 
#' @param obj a data cube proxy object (class cube)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' cat(as_json(select_bands(raster_cube(L8.col, v), c("B04", "B05"))))
#' @export
as_json <- function(obj) {
  stopifnot(is.cube(obj))
  x = libgdalcubes_cube_info(obj)
  return(x$graph)
}




#' Materialize a data cube as a netCDF file
#' 
#' This function will read chunks of a data cube and write them to a single netCDF file. The resulting
#' file uses the enhanced netCDF-4 format (for chunking and compression).
#' 
#' @seealso \url{https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf/NetCDF_002d4-Format.html#NetCDF_002d4-Format}
#' @seealso \code{\link{gdalcubes_set_ncdf_compression}} 
#' @param x a data cube proxy object (class cube)
#' @param fname output file name
#' @details The resulting netCDF file contains three dimensions (t, y, x) and bands as variables.
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' write_ncdf(select_bands(raster_cube(L8.col, v), c("B04", "B05")), fname=tempfile(fileext = ".nc"))
#' @export
write_ncdf <- function(x, fname = tempfile(pattern = "gdalcubes", fileext = ".nc")) {
  stopifnot(is.cube(x))
  
  if (.pkgenv$use_cube_cache) {
    j = as_json(x)
    if (!is.null(.pkgenv$cube_cache[[j]])
        && file.exists(.pkgenv$cube_cache[[j]])) {
      file.copy(from=.pkgenv$cube_cache[[j]], to = fname, overwrite=TRUE)
    }
    else {
      libgdalcubes_eval_cube(x, fname, .pkgenv$compression_level)
    }
  }
  else {
    libgdalcubes_eval_cube(x, fname, .pkgenv$compression_level)
  }
  invisible()
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
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
    return(libgdalcubes_dimension_values(obj, datetime_unit)) 
  }
  else if (is.cube_view(obj)) {
    if (is.null(datetime_unit)) {
      datetime_unit = ""
    }    
    return(libgdalcubes_dimension_values_from_view(obj, datetime_unit)) 
  }
  else {
    stop("obj must be either from class cube or from class cube_view")
  }
}



