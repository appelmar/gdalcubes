

serialize_function <- function(f) {
  src <- attr(f,"srcref", exact = TRUE)
  if (is.null(src))
    stop("source for given function is not available")
  return(paste(as.character(src),collapse = "\n"))
}

#' Apply a function on chunks of a data cube 
#' d
#' @details 
#' This function internally creates gdalcubes stream data cube, which streams
#' data of a chunk to stdin of a new R process. For reading data the function f typically 
#' calls \code{x <- read_chunk_as_array()} which then results in a 4 dimensional (band, time, y, x) array.
#' Similarly \code{write_stream_from_array(x)} will write a result array as a chunk in the resulting data cube.
#'
#' @param cube Source data cube
#' @param f R function to apply over all chunks
#' @return A proxy data cube object
#' @examples 
#' \dontrun{
#' f <- function() {
#'   x = read_stream_as_array()
#'   out <- reduce_time(x, function(x) {
#'     y = x[1,]
#'     mean(y, na.rm=T)})
#'   write_stream_from_array(out)
#' }
#' xstrm <- gcbs_chunk_apply(XXX, f)
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
gcbs_chunk_apply <- function(cube, f) {
  stopifnot(is.gcbs_cube(cube))
  srcfile =  tempfile(".stream_",fileext = ".R")
  cat(serialize_function(f), file = srcfile, append = FALSE)
  cmd <- paste("Rscript ", "--vanilla ", "-e ", "\"require(gdalcubes)\" ", "-e ", "\"do.call(eval(parse('", srcfile ,"')), args=list())\"", sep="")
  x = libgdalcubes_create_stream_cube(cube, cmd)
  class(x) <- c("gcbs_chunk_apply_cube", "gcbs_cube", "xptr")
  return(x)
}

is.gcbs_chunk_apply_cube  <- function(obj) {
  if(!("gcbs_chunk_apply_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


