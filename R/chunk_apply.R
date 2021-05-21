serialize_function <- function(f, args) {
  stopifnot(is.function(f))
  #src <- attr(f,"srcref", exact = TRUE)
  #if (is.null(src))
  #  stop("source for given function is not available")
  #return(paste(as.character(src),collapse = "\n"))
  if(is.null(args)){
    return(paste(deparse(f),collapse = "\n"))
  } else {
    stopifnot(is.list(args))
    f = deparse(f)
    if(grepl("=", f[1])){
      stop('Please do not use "=" in your function assignment. Use args to supply arguments insted')
    }
    for(arg in names(args)){
      if(is.character(args[[arg]])){
        val = paste0('\"', args[arg], '\"')
        f[1] = gsub(arg, paste0(arg,"=", val), f[1])
      } else {
        f[1] = gsub(arg, paste0(arg,"=", args[[arg]]), f[1])
      }
    }
    return(paste(f, collapse = "\n"))
  }
}


#' Apply an R function on chunks of a data cube
#'
#' @details
#' This function internally creates a gdalcubes stream data cube, which streams
#' data of a chunk to a new R process. For reading data, the function typically
#' calls \code{x <- read_chunk_as_array()} which then results in a 4 dimensional (band, time, y, x) array.
#' Similarly \code{write_chunk_from_array(x)} will write a result array as a chunk in the resulting data cube.
#' The chunk size of the input cube is important to control how the function will be exposed to the data cube. For example,
#' if you want to apply an R function over complete pixel time series, you must define the chunk size argument in \code{\link{raster_cube}}
#' to make sure that chunk contain the correct parts of the data.
#'
#' @param cube source data cube
#' @param f R function to apply over all chunks
#' @param args optional additional parameters to f
#' @return a proxy data cube object
#' @examples
#' \donttest{
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
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'                           srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube = select_bands(L8.cube, c("B04", "B05"))
#' f <- function() {
#'   x <- read_chunk_as_array()
#'   out <- reduce_time(x, function(x) {
#'     cor(x[1,], x[2,], use="na.or.complete", method = "kendall")
#'   })
#'   write_chunk_from_array(out)
#' }
#' L8.cor = chunk_apply(L8.cube, f)
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
chunk_apply <- function(cube, f, args=NULL) {
  stopifnot(is.cube(cube))

  funstr = serialize_function(f, args)
  funhash = libgdalcubes_simple_hash(funstr)
  srcfile =  file.path(tempdir(), paste(".stream_", funhash, ".R", sep=""))
  srcfile = gsub("\\\\", "/", srcfile) # Windows fix

  cat(funstr, file = srcfile, append = FALSE)

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


