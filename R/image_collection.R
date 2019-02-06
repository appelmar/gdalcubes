#' Load an existing image collection from a file
#'
#' This function will load an image collection from an SQLite file. These files
#' index / reference existing imagery. To create a collection from files on disk,
#' use \code{create_collection}.
#' @param path Path to an existing image collection file
#' @return An image collection proxy object, which can be used to create a data cube using \code{cube}
#' @examples 
#' L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                        ".TIF", recursive = TRUE, full.names = TRUE)
#' create_image_collection(L8_files, "L8_L1TP", out_file = "L8.db", quiet=TRUE)
#' image_collection("L8.db")
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
#' @param x image collection object
#' @param srs target spatial reference system
#' @return a list with elements \code{left}, \code{right}, \code{bottom}, \code{top}, \code{t0} (start date/time), and \code{t1} (end date/time)
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  L8.col = create_image_collection(L8_files, "L8_L1TP") 
#'  extent(L8.col,"EPSG:32618")
#'  cube_view(extent=extent(L8.col,"EPSG:32618"), 
#'            bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'            srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' @export
extent <- function(x, srs="EPSG:4326") {
  stopifnot(is.image_collection(x))
  return(libgdalcubes_image_collection_extent(x, srs))
}





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
#' This function iterates over files or GDAL dataset identifiers and extracts datetime, image_id, and band information according to a given
#' collection format.
#' 
#' @details
#' An image collection is a simple SQLite database file that indexes and references existing image files / GDAL dataset identifiers.
#' @param files A character vector with paths to image files on disk or any GDAL dataset identifiers (including virtual file systems and higher level drivers or GDAL subdatasets)
#' @param out_file Optional name of the output SQLite database file, defaults to a temporary file
#' @param format Collection format, can be either a name to use predefined formats (as output from \code{collection_formats}) or a path to a custom JSON format description file.
#' @param unroll_archives automatically convert .zip, .tar archives and .gz compressed files to GDAL virtual file system dataset identifiers (e.g. by prepending /vsizip/) and add contained files to the list of files  
#' @param quiet logical; if TRUE, do not print resulting image collection if return value is not assigned to a variable
#' @return An image collection proxy object, which can be used to create a data cube using \code{cube}
#' @examples 
#' L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"), 
#'                        ".TIF", recursive = TRUE, full.names = TRUE)
#' create_image_collection(L8_files, "L8_L1TP")
#' @export
create_image_collection <-function(files, format, out_file=tempfile(fileext = ".sqlite"), unroll_archives=TRUE, quiet=FALSE)
{
  libgdalcubes_create_image_collection(files, format, out_file, unroll_archives)
  if (quiet) {
    return(invisible(image_collection(out_file)))
  }
  return(image_collection(out_file))
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
#' @param print Should available formats and their descriptions be printed nicely?
#' @return A data.frame with columns name and description where the former describes the unique identifier that can be used in \code{create_image_collection} and the
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