#' Reduce multidimensional data over time
#' 
#' This generic function applies a reducer function over a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{reduce_time.cube}} 
#' @seealso \code{\link{reduce_time.array}} 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' reduce_time(raster_cube(L8.col, v) , "median(B02)", "median(B03)", "median(B04)")  
#'  
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' y <- reduce_time(x, function(v) {
#'   apply(v, 1, mean)
#' })
#'  
#' @export
reduce_time <- function(x, ...) {
  UseMethod("reduce_time")
}

#' Reduce multidimensional data over space
#' 
#' This generic function applies a reducer function over a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{reduce_space.cube}} 
#' @seealso \code{\link{reduce_space.array}} 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' reduce_space(raster_cube(L8.col, v) , "median(B02)")  
#' 
#' 
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' y <- reduce_space(x, function(v) {
#'   apply(v, 1, mean)
#' })
#'  
#' @export
reduce_space <- function(x, ...) {
  UseMethod("reduce_space")
}



#' Reduce a data cube over the time dimension
#' 
#' Create a proxy data cube, which applies one or more reducer functions to selected bands over pixel time series of a data cube
#'
#' @param x source data cube
#' @param expr either a single string, or a vector of strings defining which reducers will be applied over which bands of the input cube
#' @param ... optional additional expressions (if \code{expr} is not a vector)
#' @param FUN a user-defined R function applied over pixel time series (see Details)
#' @param load_pkgs logical or character; if TRUE, all currently attached packages will be attached automatically before executing FUN in spawned R processes, specific packages can alternatively be provided as a character vector.
#' @param load_env logical or environment; if TRUE, the current global environment will be restored automatically before executing FUN in spawned R processes, can be set to a custom environment.
#' @param names character vector; names of the output bands, if FUN is provided, the length of names is used as the expected number of output bands
#' @return proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.rgb.median = reduce_time(L8.rgb, "median(B02)", "median(B03)", "median(B04)")  
#' L8.rgb.median
#' 
#' \donttest{
#' plot(L8.rgb.median, rgb=3:1)
#' }
#' 
#' # user defined reducer calculating interquartile ranges
#' L8.rgb.iqr = reduce_time(L8.rgb, names=c("iqr_R", "iqr_G","iqr_B"), FUN = function(x) {
#'     c(diff(quantile(x["B04",],c(0.25,0.75), na.rm=TRUE)),
#'       diff(quantile(x["B03",],c(0.25,0.75), na.rm=TRUE)),
#'       diff(quantile(x["B02",],c(0.25,0.75), na.rm=TRUE)))
#' })
#' L8.rgb.iqr
#' \donttest{
#' plot(L8.rgb.iqr, key.pos=1)
#' }
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details 
#' 
#' The function can either apply a built-in reducer if expr is given, or apply a custom R reducer function if FUN is provided.
#' 
#' In the former case, notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments. Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median", "var", "sd", "which_min", "which_max",
#' "Q1" (1st quartile), and "Q3" (3rd quartile).
#' 
#' User-defined R reducer functions receive a two-dimensional array as input where rows correspond to the band and columns represent the time dimension. For 
#' example, one row is the time series of a specific band. FUN should always return a numeric vector with the same number of elements, which will be interpreted
#' as bands in the result cube. Notice that it is recommended to specify the names of the output bands as a character vector. If names are missing,
#' the number and names of output bands is tried to be derived automatically, which may fail in some cases. 
#' 
#' For more details and examples on how to write user-defined functions, please refer to the gdalcubes website 
#' at \url{https://gdalcubes.github.io/source/concepts/udfs.html}.
#' 
#' @export
reduce_time.cube <- function(x, expr, ..., FUN, names=NULL, load_pkgs=FALSE, load_env=FALSE) {
  stopifnot(is.cube(x))
  if (missing(expr) && missing(FUN)) {
    stop("either a expr or FUN must be provided ")
  }
  if (!missing(FUN) && !missing(expr)) {
    warning("received both expr and FUN, ignoring FUN")
    FUN = NULL
  }
  if (!missing(FUN) && !is.function(FUN)) {
    stop ("FUN must be a function")
  }
  
  
  if (!missing(expr)) {
    stopifnot(is.character(expr))
    if (length(list(...))> 0) {
      stopifnot(all(sapply(list(...), is.character)))
      expr = c(expr, unlist(list(...)))
    }
    
    # parse expr to separate reducers and bands
    reducers = gsub("\\(.*\\)", "", expr)
    bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
    stopifnot(length(reducers) == length(bands))
    x = gc_create_reduce_time_cube(x, reducers, bands, names)
    class(x) <- c("reduce_time_cube", "cube", "xptr")
    return(x)
  }
  else {
    
    if (!is.null(names)) {
      nb = length(names)
    }
    else {
      # guess number of bands from provided function
      dummy_values = matrix(rnorm(nbands(x)*10), nrow = nbands(x), ncol=10)
      rownames(dummy_values) <- names(x)
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
    cat("write_chunk_from_array(reduce_time(read_chunk_as_array(), f))", "\n", file = srcfile2, append = TRUE)
    cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
    
    x = gc_create_stream_reduce_time_cube(x, cmd, nb, names)
    class(x) <- c("reduce_time_cube", "cube", "xptr")
    return(x)
  }
  
}



#' Reduce a data cube over spatial (x,y or lat,lon) dimensions
#' 
#' Create a proxy data cube, which applies one or more reducer functions to selected bands over spatial slices of a data cube
#'
#' @param x source data cube
#' @param expr either a single string, or a vector of strings defining which reducers will be applied over which bands of the input cube
#' @param ... optional additional expressions (if \code{expr} is not a vector)
#' @param FUN a user-defined R function applied over pixel time series (see Details)
#' @param load_pkgs logical or character; if TRUE, all currently attached packages will be attached automatically before executing FUN in spawned R processes, specific packages can alternatively be provided as a character vector.
#' @param load_env logical or environment; if TRUE, the current global environment will be restored automatically before executing FUN in spawned R processes, can be set to a custom environment.
#' @param names character vector; names of the output bands, if FUN is provided, the length of names is used as the expected number of output bands
#' @return proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does).
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.b02 = select_bands(L8.cube, c("B02"))
#' L8.b02.median = reduce_space(L8.b02, "median(B02)")  
#' L8.b02.median
#' \donttest{
#' plot(L8.b02.median)
#' }
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parentheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently include "min", "max", "sum", "prod", "count", "mean", "median", "var", and "sd".
#' 
#' For more details and examples on how to write user-defined functions, please refer to the gdalcubes website 
#' at \url{https://gdalcubes.github.io/source/concepts/udfs.html}.
#' 
#' 
#' @export
reduce_space.cube <- function(x, expr, ..., FUN, names=NULL, load_pkgs=FALSE, load_env=FALSE) {
  stopifnot(is.cube(x))
  if (missing(expr) && missing(FUN)) {
    stop("either a expr or FUN must be provided ")
  }
  if (!missing(FUN) && !missing(expr)) {
    warning("received both expr and FUN, ignoring FUN")
    FUN = NULL
  }
  if (!missing(FUN) && !is.function(FUN)) {
    stop ("FUN must be a function")
  }
  
  
  if (!missing(expr)) {
    stopifnot(is.character(expr))
    if (length(list(...))> 0) {
      stopifnot(all(sapply(list(...), is.character)))
      expr = c(expr, unlist(list(...)))
    }
    
    # parse expr to separate reducers and bands
    reducers = gsub("\\(.*\\)", "", expr)
    bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
    stopifnot(length(reducers) == length(bands))
    x = gc_create_reduce_space_cube(x, reducers, bands, names)
    class(x) <- c("reduce_space_cube", "cube", "xptr")
    return(x)
  }
  else {
    if (!is.null(names)) {
      nb = length(names)
    }
    else {
      # guess number of bands from provided function
      dummy_values = matrix(rnorm(nbands(x)*10), nrow = nbands(x), ncol=10)
      rownames(dummy_values) <- names(x)
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
    cat("write_chunk_from_array(reduce_space(read_chunk_as_array(), f))", "\n", file = srcfile2, append = TRUE)
    cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
    
    x = gc_create_stream_reduce_space_cube(x, cmd, nb, names)
    class(x) <- c("reduce_space_cube", "cube", "xptr")
    return(x)
    
  }
  
}




is.reduce_time_cube  <- function(obj) {
  if(!("reduce_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.reduce_space_cube  <- function(obj) {
  if(!("reduce_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


