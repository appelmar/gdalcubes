#' Load an existing image collection from a file
#'
#' This function will load an image collection from an SQLite file. Image collection files
#' index and reference existing imagery. To create a collection from files on disk,
#' use \code{\link{create_image_collection}}.
#' @param path path to an existing image collection file
#' @return an image collection proxy object, which can be used to create a data cube using \code{\link{raster_cube}}
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
#' L8.col
#' @export
image_collection <- function(path) {
  stopifnot(file.exists(path))
  xptr <- libgdalcubes_open_image_collection(path)
  class(xptr) <- c("image_collection", "xptr")
  return(xptr)
}


is.image_collection <- function(obj) {
  if(!("image_collection" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("Image collection proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


#' Derive the spatiotemporal extent of an image collection 
#'
#' @param x image collection proxy object 
#' @param srs target spatial reference system
#' @return a list with elements \code{left}, \code{right}, \code{bottom}, \code{top}, \code{t0} (start date/time), and \code{t1} (end date/time)
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
#' extent(L8.col,"EPSG:32618")
#' cube_view(extent=extent(L8.col,"EPSG:32618"),
#'           srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' @export
extent <- function(x, srs="EPSG:4326") {
  stopifnot(is.image_collection(x))
  return(libgdalcubes_image_collection_extent(x, srs))
}




#' Print image collection information
#' 
#' Prints information about images in an image collection.
#' 
#' @param x Object of class "image_collection"
#' @param ... Further arguments passed to the generic print function
#' @param n Number of images for which details are printed 
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
#' print(L8.col)
#' @export
print.image_collection <- function(x, ..., n=6) {
  stopifnot(is.image_collection(x))
  info <- libgdalcubes_image_collection_info(x)
  cat(paste("A GDAL image collection object, referencing",nrow(info$images), "images with", nrow(info$bands), " bands\n"))
  cat("Images:\n")
  print(head(info$images, n))
  if (n < nrow(info$images)) {
    cat("[ omitted", nrow(info$images) - n , "images ]", "\n")
  }
  cat("\n")
  cat("Bands:\n")
  
  print(info$bands)
  cat("\n")
}


#' Create an image collection from a set of GDAL datasets or files
#' 
#' This function iterates over files or GDAL dataset identifiers and extracts datetime, image identifiers, and band information according to a given
#' collection format.
#' 
#' @details
#' An image collection is a simple SQLite database file that indexes and references existing image files / GDAL dataset identifiers.
#' 
#' Collections can be created in two different ways: First, if a collection format is specified (argument \code{format}), date/time, bands, and metadata are automatically
#' extracted from provided files / GDAL datasets. Second, image collections can be created without collection format by manually specifying date/time of images 
#' (argument \code{date_time}) and (optional) names of bands. In this case, however, all provided images must contain the same bands. If this is not 
#' possible for a dataset, a collection format must be used. 
#' 
#' @param files character vector with paths to image files on disk or any GDAL dataset identifiers (including virtual file systems and higher level drivers or GDAL subdatasets)
#' @param out_file optional name of the output SQLite database file, defaults to a temporary file
#' @param format collection format, can be either a name to use predefined formats (as output from \code{\link{collection_formats}}) or a path to a custom JSON format description file
#' @param unroll_archives automatically convert .zip, .tar archives and .gz compressed files to GDAL virtual file system dataset identifiers (e.g. by prepending /vsizip/) and add contained files to the list of considered files  
#' @param quiet logical; if TRUE, do not print resulting image collection if return value is not assigned to a variable
#' @param date_time vector with date/ time for files; can be of class character, Date, or POSIXct (argument is only applicable for image collections without collection format)
#' @param band_names character vector with band names, length must match the number of bands in provided files (argument is only applicable for image collections without collection format)
#' @param use_subdatasets logical; use GDAL subdatasets of provided files (argument is only applicable for image collections without collection format)
#' @return image collection proxy object, which can be used to create a data cube using \code{\link{raster_cube}}
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' @export
create_image_collection <-function(files, format=NULL, out_file=tempfile(fileext = ".sqlite"), date_time=NULL, band_names=NULL, use_subdatasets=FALSE, unroll_archives=TRUE, quiet=FALSE)
{
  
  if (is.null(format) && is.null(date_time)) {
    stop("Function expects either a collection format (argument format), or a date/time vector (argument date_time)")
  }
  if (!is.null(format) && !is.null(date_time)) {
    warning("Date/time values will be extracted using the collection format; ignoring date_time argument")
    date_time = NULL
  }
  if (is.null(date_time) && !is.null(band_names)) {
    warning("Band names will be taken from the collection format; ignoring band_names argument")
    band_names = NULL
  }
  
  if (is.null(date_time)) {
    # use collection format
    libgdalcubes_create_image_collection_from_format(files, format, out_file, unroll_archives)
  }
  else {
    # no collection format, use date_time and band_names (if given)
    if ("POSIXct" %in% class(date_time) || "Date" %in% class(date_time)) {
      date_time = format(date_time,format = "%Y-%m-%dT%H:%M:%S")
    }
    else if (!is.character(date_time)) {
      date_time = as.character(date_time)
    }
    
    if (is.null(band_names)) {
      band_names=character()
    }
      libgdalcubes_create_image_collection_from_datetime(out_file, files, date_time, use_subdatasets, band_names)
  }
  
  
  if (quiet) {
    return(invisible(image_collection(out_file)))
  }
  return(image_collection(out_file))
}

#' Add images to an existing image collection
#' 
#' This function adds provided files or GDAL dataset identifiers and to an existing image collection by extracting datetime, image identifiers, and band information according to the collection's format.
#' 
#' @param image_collection image_collection object or path to an existing collection file
#' @param files character vector with paths to image files on disk or any GDAL dataset identifiers (including virtual file systems and higher level drivers or GDAL subdatasets)
#' @param unroll_archives automatically convert .zip, .tar archives and .gz compressed files to GDAL virtual file system dataset identifiers (e.g. by prepending /vsizip/) and add contained files to the list of considered files  
#' @param out_file path to output file, an empty string (the default) will update the collection in-place, whereas images will be added to a new copy of the image collection at the given location otherwise.
#' @param quiet logical; if TRUE, do not print resulting image collection if return value is not assigned to a variable
#' @return image collection proxy object, which can be used to create a data cube using \code{\link{raster_cube}}
#' @examples 
#' L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#' L8_col = create_image_collection(L8_files[1:12], "L8_L1TP") 
#' add_images(L8_col, L8_files[13:24])
#' @export
add_images <- function(image_collection, files, unroll_archives = TRUE, out_file = "", quiet = FALSE) {
  if (is.character(image_collection)) {
    image_collection = image_collection(image_collection)
  }
  stopifnot(is.image_collection(image_collection))
  libgdalcubes_add_images(image_collection, files, unroll_archives, out_file)
  
  if (quiet) {
    return(invisible(image_collection))
  }
  return(image_collection)
}


#' List predefined image collection formats
#'
#' gdalcubes comes with some predefined collection formats e.g. to scan Sentinel 2 data. This function lists available formats  including brief descriptions.
#' 
#' @details 
#' Image collection formats define how individual files / GDAL datasets relate to an image collection, i.e., 
#' which bands they contain, to which image they belong, and how to derive aquisition date/time.
#' They are described as a set of regular expressions in a JSON file and used by gdalcubes to extract this information 
#' from the paths and/or filenames. 
#' 
#' @param print logical; should available formats and their descriptions be printed nicely, defaults to TRUE
#' @return data.frame with columns \code{name} and \code{description} where the former describes the unique identifier that can be used in \code{create_image_collection} and the
#' latter gives a brief description of the format. 
#' @examples 
#' collection_formats()
#' @export
collection_formats <-function(print=TRUE)
{
  df = libgdalcubes_list_collection_formats()
  df$name = as.character(df$name)
  df$path = as.character(df$path)
  df$description = ""
  df$tags = ""
  for (i in 1:nrow(df)) {
    x = jsonlite::read_json(df$path[i])$description
    if (!is.null(x))
      df$description[i] = x
    x = jsonlite::read_json(df$path[i])$tags
    if (!is.null(x))
      df$description[i] = paste(df$description[i], " [TAGS: ", paste(x, collapse=", "), "]", sep="") 
  }
  w = getOption("width")
  
  lw = max(nchar(df$name))
  if (lw >= w) {
    # TODO: what if current width is less than the format identifier column?
  }
  
  y = lapply(df$description, strwrap, width=w-lw-3)
  
  # TODO: header
  if (print) {
    for (i in 1:nrow(df)) {
      cat( rep(" ", lw-nchar(df$name[i])), df$name[i], " | ", y[[i]][1], sep="")
      cat("\n", sep="")
      if (length(y[[i]]) > 1) {
        for (j in 2:length(y[[i]])) {
          cat(rep(" ", lw)," | ", y[[i]][j], sep="")
          cat("\n", sep="")
        }
      }
    }
  }
  return(invisible(df[,c("name","description")]))
}



#' Download and install an image collection format from a URL
#' 
#' @param url URL pointing to the collection format JSON file
#' @param name optional name used to refer to the collection format
#' @details 
#' By default, the collection format name will be derived from the basename of the URL.
#' @examples 
#' \donttest{
#' add_collection_format(
#'    "https://raw.githubusercontent.com/appelmar/gdalcubes/dev/formats/Sentinel1_IW_GRD.json")
#' }
#' @export
add_collection_format <- function(url, name=NULL) {
  if (is.null(name)) {
    name = basename(url)
  }
  else {
    if (!endsWith(name, ".json")) {
      name = paste(name, ".json", sep="")
    }
  }
  destfile = file.path(system.file(package="gdalcubes"), "formats", name)
  download.file(url, destfile = destfile)
  if(!jsonlite::validate(readLines(destfile))) {
    file.remove(destfile)
    stop("downloaded file is not valid json")
  }
  invisible()
}
