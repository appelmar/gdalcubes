#' Load an existing image collection from a file
#'
#' This function will load an image collection from an SQLite file. These files
#' index / reference existing imagery. To create a collection from files on disk,
#' use \code{gcbs_create_collection}.
#' @param path Path to an existing image collection file
#' @return An image collection proxy object, which can be used to create a data cube using \code{gcbs_cube}
#' @examples 
#' \dontrun{gcbs_image_collection("/home/marius/github/gdalcubes/cmake-build-debug/src/test.db")}
#' @export
gcbs_image_collection <- function(path) {
  stopifnot(file.exists(path))
  xptr <- libgdalcubes_open_image_collection(path)
  class(xptr) <- c("gcbs_image_collection", "xptr")
  return(xptr)
}


is.gcbs_image_collection <- function(obj) {
  if(!("gcbs_image_collection" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("Image collection proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


#' @export
print.gcbs_image_collection <- function(x, ...) {
  stopifnot(is.gcbs_image_collection(x))
  info <- libgdalcubes_image_collection_info(x)
  cat(paste("A GDAL image collection object, referencing",nrow(info$images), "images with", nrow(info$bands), " bands\n"))
  cat("Images:\n")
  print(head(info$images))
  cat("\n")
  cat("Bands:\n")
  print(info$bands)
  cat("\n")
}


#' Create an image collection from a set of GDAL datasets or files
#' 
#' This function iterates over files or GDAL dataset identifiers and extracts datetime, image_id, and band information according to a given
#' collection format.
#' 
#' @details
#' An image collection is a simple SQLite database file that indexes and references existing image files / GDAL dataset identifiers.
#' @param files A character vector with paths to image files on disk or any GDAL dataset identifiers (including virtual file systems and higher level drivers or GDAL subdatasets)
#' @param out_file Optional name of the output SQLite database file, defaults to a temporary file
#' @param format Collection format, can be either a name to use predefined formats (as output from \code{gcbs_collection_formats}) or a path to a custom JSON format description file.
#' @param unroll_archives automatically convert .zip, .tar archives and .gz compressed files to GDAL virtual file system dataset identifiers (e.g. by prepending /vsizip/) and add contained files to the list of files  
#' @return An image collection proxy object, which can be used to create a data cube using \code{gcbs_cube}
#' @export
gcbs_create_image_collection <-function(files, format, out_file=tempfile(fileext = ".sqlite"), unroll_archives=TRUE)
{
  libgdalcubes_create_image_collection(files, format, out_file, unroll_archives)
  return(gcbs_image_collection(out_file))
}



#' List predefined image collection formats
#'
#' gdalcubes comes with some predefined collection formats e.g. to scan Sentinel 2 data. This function lists available formats and a brief description.
#' 
#' @details 
#' Image collection formats define how individul files / GDALDatasets relate to a collection, i.e., 
#' which bands they contain, to which image they belong, and how the derive date/time.
#' They are described as a set of regular expressions in a JSON file and used by gdalcubes to extract this information 
#' from the paths and/or filenames. 
#' 
#' @return A data.frame with columns name and description where the former describes the unique identifier that can be used in \code{gcbs_create_image_collection} and the
#' latter gives a brief description of the format.
#' @export
gcbs_collection_formats <-function()
{
  df = libgdalcubes_list_collection_formats()
  df$name = as.character(df$name)
  df$path = as.character(df$path)
  df$description = ""
  for (i in 1:nrow(df)) {
    x = jsonlite::read_json(df$path[i])$description
    if (!is.null(x))
      df$description[i] = x
  }
  return(df[,c("name","description")])
}