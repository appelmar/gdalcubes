#' Apply a function over (multi-band) pixels
#' 
#' This generic function applies a function on pixels of a data cube, an R array, or other classes if implemented.
#' 
#' @param x input data 
#' @param ... additional arguments passed to method implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{apply_pixel.cube}}
#' @seealso \code{\link{apply_pixel.array}} 
#' @examples 
#' 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#'               
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' apply_pixel(raster_cube(L8.col, v), "(B05-B04)/(B05+B04)", "NDVI") 
#' 
#' d <- c(4,16,128,128)
#' x <- array(rnorm(prod(d)), d)
#' y <- apply_pixel(x, function(v) {
#'   v[1] + v[2] + v[3] - v[4]
#' })
#' 
#' @export
apply_pixel <- function(x, ...) {
  UseMethod("apply_pixel")
}




#' Apply arithmetic expressions over all pixels of a data cube
#' 
#' Create a proxy data cube, which applies arithmetic expressions over all pixels of a data cube. Expressions may access band values by name.
#'
#' @param x source data cube
#' @param expr character vector with one or more arithmetic expressions (see Details)
#' @param names optional character vector with the same length as expr to specify band names for the output cube
#' @param ... not used
#' @param FUN user-defined function TODO
#' @return a proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/codeplea/tinyexpr}{tinyexpr library} to evaluate expressions in C / C++, you can look at the \href{https://github.com/codeplea/tinyexpr#functions-supported}{library documentation}
#' to see what kind of expressions you can execute. Pixel band values can be accessed by name.
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
#' L8.cube = raster_cube(L8.col, v) 
#' L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#' L8.ndvi
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
apply_pixel.cube <- function(x, expr, names=NULL, ..., FUN) {
  stopifnot(is.cube(x))
  
  
  if (missing(expr) && missing(FUN)) {
    stop("either expr or FUN must be provided ")
  }
  if (!missing(FUN) && !missing(expr)) {
    warning("received both expr and FUN, ignoring FUN")
    FUN = NULL
  }
  if (!missing(FUN) && !is.function(FUN)) {
    stop ("FUN must be a function")
  }
  
  if (!missing(expr)) {
    # apply C expression on band values
    if (is.null(names)) {
      names <- paste("band", 1:length(expr), sep="")
    }
    
    x = libgdalcubes_create_apply_pixel_cube(x, expr, names)
    class(x) <- c("apply_pixel_cube", "cube", "xptr")
    return(x)
  }
  else {
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
    srcfile1 =  tempfile(".stream_",fileext = ".R")
    srcfile1 = gsub("\\\\", "/", srcfile1) # Windows fix
    
    cat(serialize_function(FUN),  file = srcfile1, append = FALSE)
    
    srcfile2 =  tempfile(".stream_",fileext = ".R")
    srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix
    cat("require(gdalcubes)", "\n", file = srcfile2, append = FALSE)
    cat(paste("assign(\"f\", eval(parse(\"", srcfile1, "\")))", sep=""), "\n", file = srcfile2, append = TRUE)
    cat("write_chunk_from_array(apply_pixel(read_chunk_as_array(), f))", "\n", file = srcfile2, append = TRUE)
    cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
    
    x = libgdalcubes_create_stream_apply_pixel_cube(x, cmd, nb, names)
    class(x) <- c("apply_pixel_cube", "cube", "xptr")
    return(x) 

  }
  
  
}



is.apply_pixel_cube  <- function(obj) {
  if(!("apply_pixel_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




