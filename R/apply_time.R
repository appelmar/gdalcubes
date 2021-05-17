#' Apply a function over (multi-band) pixel time series
#'
#' This generic function applies a function on pixel time series of a data cube, an R array, or other classes if implemented.
#' The resulting object is expected to have the same spatial and temporal shape as the input, i.e., no reduction is performed.
#'
#' @param x input data
#' @param ... additional arguments passed to method implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{apply_time.cube}}
#' @seealso \code{\link{apply_time.array}}
#' @examples
#' # 1. input is data cube
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube = select_bands(L8.cube, c("B04", "B05"))
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI")
#'
#' # Apply a user defined R function
#' apply_time(L8.ndvi, names="NDVI_residuals",
#'    FUN=function(x) {
#'       y = x["NDVI",]
#'       if (sum(is.finite(y)) < 3) {
#'          return(rep(NA,ncol(x)))
#'       }
#'       t = 1:ncol(x)
#'       return(predict(lm(y ~ t)) -  x["NDVI",])})
#'
#' # 2. input is array
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' z <- apply_time(x, function(v) {
#'   y = matrix(NA, ncol=ncol(v), nrow=2)
#'   y[1,] = (v[1,] + v[2,]) / 2
#'   y[2,] = (v[3,] + v[4,]) / 2
#'   y
#' })
#' dim(z)
#'
#' @export
apply_time <- function(x, ...) {
  UseMethod("apply_time")
}


#' Apply a user-defined R function over (multi-band) pixel time series
#'
#' Create a proxy data cube, which applies a user-defined R function over all pixel time series of a data cube.
#' In contrast to \code{\link{reduce_time}}, the time dimension is not reduced, i.e., resulting time series
#' must have identical length as the input data cube but may contain a different number of bands / variables.
#' Example uses of this function may include time series decompositions, cumulative sums / products, smoothing, sophisticated
#' NA filling, or similar.
#'
#' @param x source data cube
#' @param names optional character vector to specify band names for the output cube
#' @param keep_bands logical; keep bands of input data cube, defaults to FALSE, i.e., original bands will be dropped
#' @param FUN user-defined R function that is applied on all pixel time series (see Details)
#' @param args optional additional arguments to FUN
#' @param ... not used
#' @return a proxy data cube object
#' @details
#' FUN receives a single (multi-band) pixel time series as a matrix with rows corresponding to bands and columns corresponding to time.
#' In general, the function must return a matrix with the same number of columns. If re result contains only a single band, it may alternatively return a vector
#' with length identical to the length of the input time series (number of columns of the input).
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube = select_bands(L8.cube, c("B04", "B05"))
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI")
#'
#' # Apply a user defined R function
#' L8.ndvi.resid = apply_time(L8.ndvi, names="NDVI_residuals",
#'    FUN=function(x) {
#'       y = x["NDVI",]
#'       if (sum(is.finite(y)) < 3) {
#'          return(rep(NA,ncol(x)))
#'       }
#'       t = 1:ncol(x)
#'       return(predict(lm(y ~ t)) -  x["NDVI",])
#'    })
#' L8.ndvi.resid
#'
#' \donttest{
#' plot(L8.ndvi.resid)
#' }
#'
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
apply_time.cube <- function(x, names=NULL, keep_bands=FALSE, FUN, args=NULL, ...) {
  stopifnot(is.cube(x))

  if (!is.function(FUN)) {
    stop ("FUN must be a function")
  }

  # apply R function on band values
  if (!is.null(names)) {
    nb = length(names)
  }
  else {
    # guess number of bands from provided function
    dummy_values = rnorm(nbands(x))
    names(dummy_values) <- names(x)
    tryCatch({
      res <- as.vector(FUN(dummy_values))
      nb <- length(res)
      # set names
      if (!is.null(names(res))) {
        names = names(res)
      }
      else {
        names = paste("band", 1:nb, sep="")
      }
    }
    , error = function(e) {
      stop("Failed to derive the length of the output from FUN automatically, please specify output band names with the correct size.")
    })
  }

  # create src file
  # TODO: load the same packages as in the current workspace? see (.packages())
  funstr = serialize_function(FUN, args)
  funhash = libgdalcubes_simple_hash(funstr)
  srcfile1 =  file.path(tempdir(), paste(".streamfun_", funhash, ".R", sep=""))
  srcfile1 = gsub("\\\\", "/", srcfile1) # Windows fix

  cat(funstr,  file = srcfile1, append = FALSE)
  srcfile2 =  file.path(tempdir(), paste(".stream_", funhash, ".R", sep=""))
  srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix

  cat("require(gdalcubes)", "\n", file = srcfile2, append = FALSE)
  cat(paste("assign(\"f\", eval(parse(\"", srcfile1, "\")))", sep=""), "\n", file = srcfile2, append = TRUE)
  cat("write_chunk_from_array(apply_time(read_chunk_as_array(), f))", "\n", file = srcfile2, append = TRUE)
  cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")

  x = libgdalcubes_create_stream_apply_time_cube(x, cmd, nb, names, keep_bands)
  class(x) <- c("apply_time_cube", "cube", "xptr")
  return(x)



}



is.apply_time_cube  <- function(obj) {
  if(!("apply_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




