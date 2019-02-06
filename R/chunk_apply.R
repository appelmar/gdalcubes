
serialize_function <- function(f) {
  src <- attr(f,"srcref", exact = TRUE)
  if (is.null(src))
    stop("source for given function is not available")
  return(paste(as.character(src),collapse = "\n"))
}

#' Apply a function on chunks of a data cube 
#' 
#' @details 
#' This function internally creates a gdalcubes stream data cube, which streams
#' data of a chunk to a new R process. For reading data the function typically 
#' calls \code{x <- read_chunk_as_array()} which then results in a 4 dimensional (band, time, y, x) array.
#' Similarly \code{write_chunk_from_array(x)} will write a result array as a chunk in the resulting data cube.
#' The chunk size of the input cube is important to control how the function will be exposed on the data cube. For example,
#' if you want to apply an R function over complete pixel time series, you must define the chunk size argument in \code{data_cube()}
#' to make sure that chunk contain the correct parts of the data. 
#' 
#'
#' @param cube Source data cube
#' @param f R function to apply over all chunks
#' @return A proxy data cube object
#' @examples 
#' \dontrun{
#' f <- function() {
#'   x = read_chunk_as_array()
#'   out <- reduce_time(x, function(x) {
#'     y = x[1,]
#'     mean(y, na.rm=T)})
#'   write_chunk_from_array(out)
#' }
#' xstrm <- chunk_apply(XXX, f)
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
chunk_apply <- function(cube, f) {
  stopifnot(is.cube(cube))
  srcfile =  tempfile(".stream_",fileext = ".R")
  srcfile = gsub("\\\\", "/", srcfile) # Windows fix
  cat(serialize_function(f), file = srcfile, append = FALSE)
  
  cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", "-e ", "\"require(gdalcubes)\" ", "-e ", "\"do.call(eval(parse('", srcfile ,"')), args=list())\"", sep="")
  x = libgdalcubes_create_stream_cube(cube, cmd)
  class(x) <- c("chunk_apply_cube", "cube", "xptr")
  return(x)
}

is.chunk_apply_cube  <- function(obj) {
  if(!("chunk_apply_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


