#' Read chunk data of a data cube from stdin or a file
#' 
#' This function can be used within function passed to \code{\link{chunk_apply}} in order to read a data cube chunk as a four-dimensional R array.
#' It works only for R processes, which have been started from the gdalcubes C++ library. 
#' The resulting array has dimensions band, time, y, x (in this order).
#'
#' @note Call this function ONLY from a function passed to \code{\link{chunk_apply}}.
#' @note This function only works in R sessions started from gdalcubes streaming.
#'
#' @param with.dimnames if TRUE, the resulting array will contain dimnames with coordinates, datetime, and band names
#' 
#' @return four-dimensional array
#' 
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
#' plot(L8.cor, zlim=c(0,1), key.pos=1)
#' }
#' @export
read_chunk_as_array <-function(with.dimnames=TRUE) {
  if(!.is_streaming()) {
    stop("This function only works in streaming mode")
  }
  
  if (Sys.getenv("GDALCUBES_STREAMING_FILE_IN") != "") {
    f <- file(Sys.getenv("GDALCUBES_STREAMING_FILE_IN"), "rb")
  }
  else {
    f <-file("stdin", "rb")
  }
 
  on.exit(close(f))
  s <- readBin(f, integer(), n=4)
  if (prod(s) == 0) {
    warning("gdalcubes::read_stream_as_array(): received empty chunk.")
    return(NULL)
  }
  bandnames <- character(s[1])
 
  for (i in 1:s[1]) {
    nchars= readBin(f, integer(), n=1)
    bandnames[i] = readChar(f, nchars = nchars)
  }
  
  dims <- readBin(f, double(), n=sum(s[2:4]))
  proj.length = readBin(f, integer(), n=1)
  proj = readChar(f, nchars = proj.length)
  buf <- readBin(f, double(), n = prod(s))
  #message("Input chunk size is ", s[1], "x",s[2], "x",s[3], "x",s[4])
  #message(paste("RECEIVED", length(buf), "/", prod(s) , "values"))

  # row major -> column major
  x <- array(buf, dim=s)
  dim(x) <- rev(dim(x))
  if (with.dimnames) {
    dnames <- list(band=bandnames,
                   datetime=dims[1:s[2]],
                   y = dims[(s[2]+1):(s[2]+s[3])],
                   x = dims[(s[2]+s[3]+1):(s[2]+s[3]+s[4])])
    dimnames(x) <- rev(dnames)
  }
  return(aperm(x,4:1))
}



#' Write chunk data of a cube to stdout or a file
#' 
#' 
#' This function can be used within function passed to \code{\link{chunk_apply}} in order to pass four-dimensional R arrays as a
#' data cube chunk to the gdalcubes C++ library. It works only for R processes, which have been started from the gdalcubes C++ library. 
#' The input array must have dimensions band, time, y, x (in this order).
#' 
#' @note Call this function ONLY from a function passed to \code{\link{chunk_apply}}.
#' @note This function only works in R sessions started from gdalcubes streaming.
#'
#' @param v four-dimensional array with dimensions band, time, y, and x
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
#' plot(L8.cor, zlim=c(0,1), key.pos=1)
#' }
#' @export
write_chunk_from_array <- function(v) {
  if(!.is_streaming()) {
    stop("This function only works in streaming mode")
  }
  v = aperm(v,c(4,3,2,1))
  dim(v) <- rev(dim(v))
  stopifnot(length(dim(v)) == 4)
  
  if (Sys.getenv("GDALCUBES_STREAMING_FILE_OUT") != "") {
    f <- file(Sys.getenv("GDALCUBES_STREAMING_FILE_OUT"), "wb")
  }
  else { # this does not work on Windows, C++ part makes sure that $GDALCUBES_STREAMING_FILE_OUT is set for Windows  
    f <- pipe("cat", "wb")
  }
  on.exit(close(f))
  s <- dim(v) 
  writeBin(as.integer(s), f)
  writeBin(as.double(v), f)
  flush(f)
}



#' Apply a function over time and bands in a four-dimensional (band, time, y, x) array and reduce time dimension
#' 
#' @param x four-dimensional input array with dimensions band, time, y, x (in this order)
#' @param FUN function which receives one time series in a two-dimensional array with dimensions bands, time as input
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result.
#' @examples
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' # reduce individual bands over pixel time series
#' y <- reduce_time(x, function(v) {
#'   apply(v, 1, mean)
#' })
#' dim(y)
#' @note This is a helper function that uses the same dimension ordering as gdalcubes streaming. It can be used to simplify 
#' the application of R functions e.g. over time series in a data cube.
#' @export
reduce_time.array <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(3,4), FUN, ...)
  if (length(dim(res)) == 2) {
    dim(res) <- c(1,1,dim(res)[1],dim(res)[2])
  }
  else if (length(dim(res)) == 3) {
    dim(res) <- c(dim(res)[1],1,dim(res)[2],dim(res)[3]) # if a vector of length n is returned, elements are interpreted as new bands of the output
  }
  return(res)
}


#' Apply a function over pixels in a four-dimensional (band, time, y, x) array
#' 
#' @param x four-dimensional input array with dimensions band, time, y, x (in this order)
#' @param FUN function that receives a vector of band values in a one-dimensional array
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result.
#' @examples
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' y <- apply_pixel(x, function(v) {
#'   v[1] + v[2] + v[3] - v[4]
#' })
#' dim(y)
#' @note This is a helper function that uses the same dimension ordering as gdalcubes. It can be used to simplify 
#' the application of R functions e.g. over time series in a data cube.
#' @export
apply_pixel.array <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(2,3,4), FUN, ...)
  dim(res) <- c(rep(1, 4-length(dim(res))), dim(res))
  return(res)
}


#' Apply a function over space and bands in a four-dimensional (band, time, y, x) array and reduce
#' spatial dimensions
#' 
#' @param x four-dimensional input array with dimensions band, time, y, x (in this order)
#' @param FUN function which receives one spatial slice in a three-dimensional array with dimensions bands, y, x as input
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result.
#' @note This is a helper function that uses the same dimension ordering as gdalcubes streaming. It can be used to simplify 
#' the application of R functions e.g. over spatial slices in a data cube.
#' @examples
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' # reduce individual bands over spatial slices 
#' y <- reduce_space(x, function(v) {
#'   apply(v, 1, mean)
#' })
#' dim(y)
#' @export
reduce_space.array <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(2), FUN, ...)
  if (is.null(dim(res)) || length(dim(res)) == 1) {
    dim(res) <- c(length(res),1,1,1) # if a vector of length n is returned, elements are interpreted as new bands of the output
    # FIXME: automatically find out whether values refer to band or time
  }
  else {
    dim(res) <- c(dim(res)[1],dim(res)[2],1,1) # if a vector of length n is returned, elements are interpreted as new bands of the output
  }
  return(res)
}









#' Apply a function over pixel time series in a four-dimensional (band, time, y, x) array
#' 
#' @param x four-dimensional input array with dimensions band, time, y, x (in this order)
#' @param FUN function that receives a vector of band values in a one-dimensional array
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a matrix (or vector if result has only one band) where rows are interpreted as new bands and columns represent time.
#' @examples
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' z <- apply_time(x, function(v) {
#'   y = matrix(NA, ncol=ncol(v), nrow=2)
#'   y[1,] = (v[1,] + v[2,]) / 2
#'   y[2,] = (v[3,] + v[4,]) / 2
#'   y
#' })
#' dim(z)
#' @note This is a helper function that uses the same dimension ordering as gdalcubes. It can be used to simplify 
#' the application of R functions e.g. over time series in a data cube.
#' @export
apply_time.array <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(3,4), FUN = FUN, ...)
  dim(res) <- c(dim(res)[1] / dim(x)[2], dim(x)[2], dim(x)[3], dim(x)[4])
  return(res)
}












