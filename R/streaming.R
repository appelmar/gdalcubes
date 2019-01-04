#' Read chunk data of a cube from stdin
#' 
#' gdalcubes stream cubes write chunk data to stdin of external processes (such as R)
#' This function can be used to read the data from stdin as four-dimensional R array with dimensions
#' (band, time, y, x).
#' 
#' @note This function only works in R sessions started from gdalcubes streaming
#' 
#' @param with.dimnames if TRUE, the resulting array will contain dimnames with coordinates, datetime, and band names
#' @return four-dimensional array
#' @export
gcbs_read_stream_as_array <-function(with.dimnames=TRUE) {
  if(!.is_streaming()) {
    stop("This function only works in streaming mode")
  }
  f <-file('stdin', 'rb')
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



#' Write chunk data of a cube to stdout
#' 
#' gdalcubes stream cubes write chunk data to stdin of external processes (such as R) and 
#' reads result data from stdout. 
#' This function can be used to write an R array as a reuslt chunk to stdout.
#' 
#' @note This function only works in R sessions started from gdalcubes streaming
#' 
#' @param v A four-dimensional array with dimensions band, time, y, and x
#' @export
gcbs_write_stream_from_array <- function(v) {
  if(!.is_streaming()) {
    stop("This function only works in streaming mode")
  }
  v = aperm(v,c(4,3,2,1))
  dim(v) <- rev(dim(v))
  stopifnot(length(dim(v)) == 4)
  f <- pipe("cat", "wb")
  on.exit(close(f))
  s <- dim(v) # test
  writeBin(as.integer(s), f)
  writeBin(as.double(v), f)
  flush(f)
}



#' Apply a function over time and bands in a four-dimensional (band, time, y, x) array
#' 
#' @param x four-dimensional input bands with dimension order band, time, y, x
#' @param FUN function which receives one time series in a two-dimensional array with dimensions bands, time as input
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result
#' @examples
#' load(system.file("extdata","sample_chunk.Rdata", package="gdalcubes"))
#' reduce_time(sample_chunk, function(x) {
#'   ndvi <- (x[8,]-x[4,])/(x[8,]+x[4,])
#'   return(c(min(ndvi, na.rm=TRUE),max(ndvi, na.rm=T)))
#' })
#' @note This is a helper function that uses the same dimension ordering as gdalcubes streaming. It can be used to simplify 
#' the application of R functions e.g. over time series in a data cube.
#' @export
reduce_time <- function(x, FUN, ...) {
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
#' @param x four-dimensional input bands with dimension order band, time, y, x
#' @param FUN function which receives a vector of band values in a one-dimensional array
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result
#' @examples
#' load(system.file("extdata","sample_chunk.Rdata", package="gdalcubes"))
#' y = apply_pixel(sample_chunk, function(x) {
#'  ndvi <- (x[8]-x[4])/(x[8]+x[4])
#'  return(c(ndvi=(x[8]-x[4])/(x[8]+x[4]), nir=x[8]))
#' })
#' @note This is a helper function that uses the same dimension ordering as gdalcubes streaming. It can be used to simplify 
#' the application of R functions e.g. over time series in a data cube.
#' @export
apply_pixel <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(2,3,4), FUN, ...)
  dim(res) <- c(rep(1, 4-length(dim(res))), dim(res))
  return(res)
}


#' Apply a function over space and bands in a four-dimensional (band, time, y, x) array
#' 
#' @param x four-dimensional input bands with dimension order band, time, y, x
#' @param FUN function which receives one spatial slice in a three-dimensional array with dimensions bands, y, x as input
#' @param ... further arguments passed to FUN
#' @details 
#' FUN is expected to produce a numeric vector (or scalar) where elements are interpreted as new bands in the result
#' @note This is a helper function that uses the same dimension ordering as gdalcubes streaming. It can be used to simplify 
#' the application of R functions e.g. over spatial slices in a data cube.
#' @examples
#' load(system.file("extdata","sample_chunk.Rdata", package="gdalcubes"))
#' y = reduce_space(sample_chunk, function(x) {
#'  ndvi <- (x[8,,]-x[4,,])/(x[8,,]+x[4,,])
#'  return(c(min(ndvi, na.rm=TRUE),max(ndvi, na.rm=T)))
#' })
#' @export
reduce_space <- function(x, FUN, ...) {
  stopifnot(is.array(x))
  stopifnot(length(dim(x))==4)
  res <- apply(x, c(2), FUN, ...)
  dim(res) <- c(dim(res)[1],dim(res)[2],1,1) # if a vector of length n is returned, elements are interpreted as new bands of the output
  return(res)
}



