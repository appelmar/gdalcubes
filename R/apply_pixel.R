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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#'               
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' apply_pixel(raster_cube(L8.col, v), "(B05-B04)/(B05+B04)", "NDVI") 
#' 
#' 
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
#' @param keep_bands logical; keep bands of input data cube, defaults to FALSE, i.e. original bands will be dropped
#' @param ... not used
#' @param FUN user-defined R function that is applied on all pixels (see Details)
#' @param load_pkgs logical or character; if TRUE, all currently attached packages will be attached automatically before executing FUN in spawned R processes, specific packages can alternatively be provided as a character vector.
#' @param load_env logical or environment; if TRUE, the current global environment will be restored automatically before executing FUN in spawned R processes, can be set to a custom environment.
#' @return a proxy data cube object
#' @details 
#' 
#' The function can either apply simple arithmetic C expressions given as a character vector (expr argument), or apply a custom R reducer function if FUN is provided.
#' 
#' In the former case, gdalcubes uses the \href{https://github.com/codeplea/tinyexpr}{tinyexpr library} to evaluate expressions in C / C++, you can look at the \href{https://github.com/codeplea/tinyexpr#functions-supported}{library documentation}
#' to see what kind of expressions you can execute. Pixel band values can be accessed by name. Predefined variables that can be used within the expression include integer pixel indexes (\code{ix}, \code{iy}, \code{it}), and 
#' pixel coordinates (\code{left}, \code{right}, \code{top}, \code{bottom}), \code{t0}, \code{t1}), where the last two values are provided seconds since epoch time.
#'  
#' FUN receives values of the bands from one pixel as a (named) vector and should return a numeric vector with identical length for all pixels. Elements of the
#' result vectors will be interpreted as bands in the result data cube. Notice that by default, since FUN is executed in a separate
#' R process, it cannot access any variables from outside and required packages must be loaded within FUN. To restore the current environment and
#' automatically load packages, set \code{load_env} and/or \code{load_pkgs} to \code{TRUE}.
#' 
#' For more details and examples on how to write user-defined functions, please refer to the gdalcubes website 
#' at \url{https://gdalcubes.github.io/source/concepts/udfs.html}.
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' # 1. Apply a C expression
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#' L8.ndvi
#' 
#' \donttest{
#' plot(L8.ndvi)
#' }
#' 
#' # 2. Apply a user defined R function
#' L8.ndvi.noisy = apply_pixel(L8.cube, names="NDVI_noisy", 
#'    FUN=function(x) {
#'        rnorm(1, 0, 0.1) + (x["B05"]-x["B04"])/(x["B05"]+x["B04"])
#'    })
#' L8.ndvi.noisy
#' 
#' 
#' \donttest{
#' plot(L8.ndvi.noisy)
#' }
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
apply_pixel.cube <- function(x, expr, names=NULL, keep_bands=FALSE, ..., FUN, load_pkgs=FALSE, load_env=FALSE) {
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
    
    x = gc_create_apply_pixel_cube(x, expr, names, keep_bands)
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
    
    
    if (is.logical(load_env)) {
      if (load_env) {
        load_env = .GlobalEnv
      }
      else
        load_env = NULL
    }
    if (!is.null(load_env)) {
      if (!is.environment(load_env)) {
        warning("Expected either FALSE/TRUE or environment for load_env; parameter will be set to FALSE.")
        load_env = NULL
      }
    }
    
    if (is.logical(load_pkgs)) {
      if (load_pkgs) {
        load_pkgs = .packages()
      }
      else {
        load_pkgs = NULL
      }
    }
    if (!is.null(load_pkgs)) {
      if (!is.character(load_pkgs)) {
        warning("Expected either FALSE/TRUE or character vector for load_pkgs; parameter will be set to FALSE.")
        load_pkgs = NULL
      }
    }
    
    
    # create src file
    # TODO: load the same packages as in the current workspace? see (.packages())
    funstr = serialize_function(FUN)
    funhash = gc_simple_hash(funstr)
    srcfile1 =  file.path(tempdir(), paste(".streamfun_", funhash, ".R", sep=""))
    srcfile1 = gsub("\\\\", "/", srcfile1) # Windows fix
    
    cat(funstr,  file = srcfile1, append = FALSE)
    srcfile2 =  file.path(tempdir(), paste(".stream_", funhash, ".R", sep=""))
    srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix

    # support custom library paths
    cat(paste0(".libPaths(",  paste(deparse(.libPaths()),collapse=""), ")\n"), file = srcfile2, append = FALSE) 
    
    cat("require(gdalcubes)", "\n", file = srcfile2, append = TRUE)
    if (!is.null(load_pkgs)) {
      cat(paste0("require(", load_pkgs,")",collapse  = "\n"), "\n", file = srcfile2, append = TRUE) 
    }
    if (!is.null(load_env)) {
      if (sum(sapply(ls(envir = load_env), FUN = function(x) {object.size(get(x, envir = load_env))})) > 100*1024^2) {
        warning("The current environment seems to be rather large (> 100 Mb), if this results in reduced performance, please consider setting load_env = FALSE.")
      }
      envfile = tempfile(pattern="renv_", fileext = ".rda")
      save(list = ls(envir = load_env),file = envfile, envir = load_env)
      cat(paste0("load(\"", envfile, "\")"), "\n", file = srcfile2, append = TRUE)
    }
    cat(paste("assign(\"f\", eval(parse(\"", srcfile1, "\")))", sep=""), "\n", file = srcfile2, append = TRUE)
    cat("write_chunk_from_array(apply_pixel(read_chunk_as_array(), f))", "\n", file = srcfile2, append = TRUE)
    cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
    
    x = gc_create_stream_apply_pixel_cube(x, cmd, nb, names, keep_bands)
    class(x) <- c("apply_pixel_cube", "cube", "xptr")
    return(x) 

  }
  
  
}



is.apply_pixel_cube  <- function(obj) {
  if(!("apply_pixel_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




