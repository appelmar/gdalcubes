#' Select a data cube band by name
#' 
#' Select a data cube band by name
#' @name gdalcubes_selection
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P3M", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' L8.red = L8.cube$B04
#' 
#' \donttest{
#' plot(L8.red)
#' }
#' @param x source data cube
#' @param name character; name of selected band
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
"$.cube" = function(x, name) {
  stopifnot(is.cube(x))
  return(select_bands(x, name))
}


#' Extract a subset of a data cube
#' 
#' Extract a subset of a data cube
#' @name gdalcubes_selection
#' @details 
#' The \code{[]} operator allows for flexible subsetting of data cubes by date, datetime,  
#' bounding box, spatial points, and band names. Depending on the arguments, it supports slicing (
#' selecting one element of a dimension) and cropping (selecting a subinterval of a dimension) and combinations
#' thereof (e.g., selecting a spatial window and a temporal slice). Dimension subsets can 
#' be specified by integer indexes or coordinates / datetime values. Arguments are matched by type and order.
#' For example, if the first argument is a length-two vector of type Date, the function will realize that this 
#' is for subsetting the time dimension. If needed, arguments are treated in the order band, time, y, x. 
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"))
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931, t0="2018-01-01", t1="2018-12-31"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1D", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' 
#' L8.cube[c("B05","B04")] # select bands
#' L8.cube[as.Date(c("2018-01-10", "2018-01-20"))] # crop by time
#' L8.cube[as.Date("2018-01-10")] # slice by time
#' L8.cube["B05", "2018-01-10"] # select bands and slice by time
#' L8.cube["B05", c("2018-01-10","2018-01-17")] # select bands and crop by time
#' L8.cube[, c("2018-01-10","2018-01-17")] # crop by time
#'
#' # crop by space (coordinates and integer indexes respectively)
#' L8.cube[list(left=388941.2 + 1e5, right=766552.4 - 1e5, bottom=4345299 + 1e5, top=4744931 - 1e5)]
#' L8.cube[,,c(1,100), c(1,100)] 
#' 
#' L8.cube[,c(1,2),,] # crop by time (integer indexes)
#' 
#' # select by spatial point or bounding box
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   s = sf::st_sfc(sf::st_point(c(500000, 4500000)), crs = "EPSG:32618")
#'   L8.cube[s]
#' 
#'   bbox =  sf::st_bbox(c(xmin = 388941.2 + 1e5, xmax = 766552.4 - 1e5,
#'                    ymax = 4744931 - 1e5, ymin = 4345299 + 1e5), crs = sf::st_crs(32618))
#'   L8.cube[bbox]
#' }
#' 
#' @param cube source data cube
#' @param ib first selector (optional), object of type character, list, Date, POSIXt, numeric, \link[sf]{st_bbox}), or \link[sf]{st_sfc}), see Details and examples
#' @param it second selector (optional), see \code{ib}
#' @param iy third selector (optional), see \code{ib}
#' @param ix fourth selector (optional), see \code{ib}
#' @param ... further arguments, not used
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
"[.cube" = function(cube, ib=TRUE, it=TRUE, iy=TRUE, ix=TRUE, ...) {
  stopifnot(is.cube(cube))
  
  args <- list(...)
  args = c(list(ib), list(it), list(iy), list(ix), args)
  
  has_bands = NULL
  has_time  = NULL
  has_y     = NULL
  has_x     = NULL
  
  x_is_int  = FALSE
  y_is_int  = FALSE
  t_is_int  = FALSE
  
  
  subset_bands <- function(X) {
    if (!is.null(has_bands)) {
      stop("Subset of band / variable dimension has already been specified")
    }
    if (is.numeric(X)) {
      if (all(X %% 1 == 0)) {
        X = names(cube)[X]
      }
    }
    if (!is.character(X)) {
      stop("Variables / bands of a data cube must be defined by integer index or their name")
    }
    has_bands = X
    return(has_bands)
  }
  subset_time <- function(X) {
    if (!is.null(has_time)) {
      stop("Subset of time dimension has already been specified")
    }
    nt = length(X)
    if (nt > 2) {
      warning(paste("Time dimension subset defined by more than two values, using range() to extract minimum and maximum as crop limits"))
      X = range(X)
    }
    if (inherits(X, "Date")) {
      has_time = format(X, "%Y-%m-%d")
    }
    else if (inherits(X, "POSIXt")) {
      has_time = format(X, "%Y-%m-%dT%H:%M:%S")
    }
    else if (is.numeric(X)) {
      if (all(X %% 1 == 0)) {
        has_time = X 
        t_is_int = TRUE
      }
      else  {
        stop("Specified time dimension subset has invalid type")
      }
    }
    else if (is.character(X)){
      has_time = X # e.g. for character YYYY-MM format´
    }
    else {
      stop("Specified time dimension subset has invalid type")
    }
    return(list(has_time, t_is_int))
  }
  subset_space <- function(X) {
    if (!is.null(has_x) || !is.null(has_y)) {
      stop("Subset of x and/or y dimension has already been specified")
    }
    if (inherits(X, "bbox")) {
      if (!requireNamespace("sf", quietly = TRUE)) {
        stop("package sf required for subsetting by spatial point; please install sf first")
      }
      sf::st_bbox(sf::st_transform(sf::st_as_sfc(X), srs(cube)))
      has_x = c( X["xmin"], X["xmax"])
      has_y = c( X["ymin"], X["ymax"])
    }
    else if (inherits(X, "sfc")) {
      if (!requireNamespace("sf", quietly = TRUE)) {
        stop("package sf required for subsetting by spatial point; please install sf first")
      }
      if (length(X) == 1) {
        if (sf::st_is(X, "POINT")) {
          X = sf::st_transform(X, srs(cube))
          has_x = sf::st_coordinates(X)[1]
          has_y = sf::st_coordinates(X)[2]
        }
        else {
          stop("Spatial selection by [] only supports bbox and point features")
        }
      }
      else {
        stop("Expected a single spatial feature")
      }
    }
    else if (is.list(X)) {
      if (!is.null(X$left)) {
        has_x = c(X$left, NA)
      }
      if (!is.null(X$right)) {
        if (is.null(has_x)) {
          has_x = c(NA, X$right)
        }
        else {
          has_x[2] = X$right
        }
      }
      if (!is.null(X$bottom)) {
        has_y = c(X$bottom, NA)
      }
      if (!is.null(X$top)) {
        if (is.null(has_y)) {
          has_y = c(NA, X$top)
        }
        else {
          has_y[2] = X$top
        }
      }
      if (!is.null(has_y) && is.null(has_x)) {
        has_x = c(NA, NA)
      }
      else if (is.null(has_y) && !is.null(has_x)) {
        has_y = c(NA, NA)
      }
      if (!is.null(X$x) && !is.null(X$y)) {
        if (is.null(has_x) && is.null(has_y)) {
          has_x = X$x
          has_y = X$x
        }
      }
      if (is.null(has_y) && is.null(has_x)) {
        warning("Definition of spatial subset is incomplete / invalid; ignoring subsetting by space")
      }
    }
    return(list(has_x, has_y, x_is_int, y_is_int))
  }
  subset_x <- function(X) {
    if (!is.null(has_x)) {
      stop("Subset of x dimension has already been specified")
    }
    if (!is.numeric(X)) {
      stop("Definition of x dimension subset expects a numeric vector")
    }
    if (length(X) > 2) {
      warning(paste("x dimension subset defined by more than two values, using range() to extract minimum and maximum as crop limits"))
      X = range(X)
    }
    if (all((X %% 1 == 0) & (X >= 1) & (X < dim(cube)[3]))) {
      # assume integer indexes
      x_is_int = TRUE
    }
    has_x = X 
    return(list(has_x, x_is_int))
  }
  subset_y <- function(X) {
    if (!is.null(has_y)) {
      stop("Subset of y dimension has already been specified")
    }
    if (!is.numeric(X)) {
      stop("Definition of y dimension subset expects a numeric vector")
    }
    if (length(X) > 2) {
      warning(paste("y dimension subset defined by more than two values, using range() to extract minimum and maximum as crop limits"))
      X = range(X)
    }
    if (all((X %% 1 == 0) & (X >= 1) & (X < dim(cube)[2]))) {
      # assume integer indexes
      y_is_int = TRUE
    }
    has_y = X 
    return(list(has_y, y_is_int))
  }
  
  
  
  
  
  
  
  
  # first, use dimensions specified by name
  arg_names = names(args) 
  if (is.null(arg_names)) {
    arg_names = rep("",length(args))
  }
  if (length(which(arg_names  != "")) > 0) {
    for (i in which(arg_names != "")) {
      name = names(args)[i] 
      if (name %in% c("t", "time", "datetime")) {
        a = subset_time(args[[i]])
        has_time = a[[1]]
        t_is_int = a[[2]]
      }
      else if (name %in% c("b", "var", "band", "variable")) {
        has_bands = subset_bands(args[[i]])
      }
      else if (name %in% c("s", "bbox", "extent", "xy", "space")) {
        a = subset_space(args[[i]])
        has_x = a[[1]]
        has_y = a[[2]]
        x_is_int = a[[3]]
        y_is_int = a[[4]]
      }
      else if (name %in% c("x", "lon")) {
        a = subset_x(args[[i]])
        has_x = a[[1]]
        x_is_int = a[[2]]
      }
      else if (name %in% c("y", "lat")) {
        a = subset_y(args[[i]])
        has_y = a[[1]]
        y_is_int = a[[2]]
      }
      else {
        stop(paste("Unknown dimension name '", name ,"'", sep=""))
      } 
    }
  }
  
  # now guess dimension for other arguments based on type and what is left
  count_missing_args = 0 # e.g. skipped dimensions as in x[,c("2018-01-01", "2019-01-01")]
  if (length(which(arg_names  == "")) > 0) {
    for (i in which(arg_names == "")) {
      if (is.logical(args[[i]]) && args[[i]] == TRUE) { # default values
        count_missing_args = count_missing_args + 1
        next
      }
      if (inherits(args[[i]], "bbox")) {
        a = subset_space(args[[i]])
        has_x = a[[1]]
        has_y = a[[2]]
        x_is_int = a[[3]]
        y_is_int = a[[4]]
      }
      else if (inherits(args[[i]], "sfc")) {
        a = subset_space(args[[i]])
        has_x = a[[1]]
        has_y = a[[2]]
        x_is_int = a[[3]]
        y_is_int = a[[4]]
      }
      else if (is.list(args[[i]])) {
        a = subset_space(args[[i]])
        has_x = a[[1]]
        has_y = a[[2]]
        x_is_int = a[[3]]
        y_is_int = a[[4]]
      }
      else if (inherits(args[[i]], "Date")) {
        a = subset_time(args[[i]])
        has_time = a[[1]]
        t_is_int = a[[2]]
      }
      else if (inherits(args[[i]], "POSIXt")) {
        a = subset_time(args[[i]])
        has_time = a[[1]]
        t_is_int = a[[2]]
      }
      else if (is.character(args[[i]])) {
        if (is.null(has_bands) && count_missing_args == 0) {
          has_bands = subset_bands(args[[i]])
        }
        else if (is.null(has_time)) {
          a = subset_time(args[[i]])
          has_time = a[[1]]
          t_is_int = a[[2]]
        }
        else {
         stop("Unknown dimension to subset by character vector, band and time dimensions have already been defined")
        }
      }
      else if (is.numeric(args[[i]])) {
        if (is.null(has_bands) && count_missing_args == 0) {
          has_bands = subset_bands(args[[i]])
        }
        else if (is.null(has_time) && count_missing_args <= 1) {
          a = subset_time(args[[i]])
          has_time = a[[1]]
          t_is_int = a[[2]]
        }
        else if (is.null(has_y) && count_missing_args <= 2) {
          a = subset_y(args[[i]])
          has_y = a[[1]]
          y_is_int = a[[2]]
        }
        else if (is.null(has_x) && count_missing_args <= 3) {
          a = subset_x(args[[i]])
          has_x = a[[1]]
          x_is_int = a[[2]]
        }
        else {
          stop("Unknown dimension to subset by numeric vector")
        }
      }
    }
  }
  
  # now, derive corresponding gdalcubes operation(s)
  out = cube
  if (!is.null(has_bands)) {
    out = .copy_cube(out)
    out = select_bands(out, has_bands)
  }
  
  if (length(has_x) != length(has_x)) {
    stop("spatial subset is invalid, make sure to define both X, and Y, or jointly as a bounding box, spatial point, or as a list")
  }
  if (x_is_int != y_is_int) {
    stop("spatial subset is invalid, make sure to define both X, and Y either as integer indexes or numeric coordinates")
  }
  
  # slicing
  if (length(has_x) == 1) {
    # spatial slicing
    if (x_is_int) {
      out = slice_space(out, i = c(has_x, has_y))
    }
    else {
      out = slice_space(out, loc = c(has_x, has_y))
    }
  }
  if (length(has_time) == 1) {
    # spatial slicing
    if (t_is_int) {
      out = slice_time(out, it = has_time)
    }
    else {
      out = slice_time(out, datetime = has_time)
    }
  }
  
  if (length(has_x) == 2 && length(has_time) == 2) {
    # spatiotemporal cropping
    if (!x_is_int && !t_is_int) {
      ee = list()
      if (is.finite(has_x[1])) {
        ee$left = has_x[1]
      }
      if (is.finite(has_x[2])) {
        ee$right = has_x[2]
      }
      if (is.finite(has_y[1])) {
        ee$bottom = has_y[1]
      }
      if (is.finite(has_y[2])) {
        ee$top= has_y[2]
      }
      if (is.finite(has_time[1])) {
        ee$t0 = has_time[1]
      }
      if (is.finite(has_time[2])) {
        ee$t1 = has_time[2]
      }
      out = crop(out, extent = ee)
    }
    else if (x_is_int && t_is_int) {
      out = crop(out, iextent = list(x = c(has_x[1], has_x[2]), y = c(has_y[1], has_y[2]),
                               t = c(has_time[1], has_time[2]) ))
    }
    else if (x_is_int && !t_is_int) {
      ee = list()
      if (is.finite(has_time[1])) {
        ee$t0 = has_time[1]
      }
      if (is.finite(has_time[2])) {
        ee$t1 = has_time[2]
      }
      out = crop(out, extent = ee)
      out = crop(out, iextent = list(x = c(has_x[1], has_x[2]), y = c(has_y[1], has_y[2]))) 
    }
    else if (!x_is_int && t_is_int) {
      ee = list()
      if (is.finite(has_x[1])) {
        ee$left = has_x[1]
      }
      if (is.finite(has_x[2])) {
        ee$right = has_x[2]
      }
      if (is.finite(has_y[1])) {
        ee$bottom = has_y[1]
      }
      if (is.finite(has_y[2])) {
        ee$top= has_y[2]
      }
      out = crop(out, extent = ee)
      out = crop(out, iextent = list(t = c(has_time[1], has_time[2])))
    }
    
  }
  else if (length(has_x) == 2) {
    # spatial cropping only
    if (x_is_int) {
      out = crop(out, iextent = list(x = c(has_x[1], has_x[2]), y = c(has_y[1], has_y[2]))) 
    }
    else {
      ee = list()
      if (is.finite(has_x[1])) {
        ee$left = has_x[1]
      }
      if (is.finite(has_x[2])) {
        ee$right = has_x[2]
      }
      if (is.finite(has_y[1])) {
        ee$bottom = has_y[1]
      }
      if (is.finite(has_y[2])) {
        ee$top= has_y[2]
      }
      out = crop(out, extent = ee)
    }
  }
  else if (length(has_time) == 2) {
    # temporal cropping only
    if (t_is_int) {
      out = crop(out, iextent = list(t = c(has_time[1], has_time[2])))
    }
    else {
      ee = list()
      if (!is.na(has_time[1]) && has_time[1] != "") {
        ee$t0 = has_time[1]
      }
      if (!is.na(has_time[2]) && has_time[2] != "") {
        ee$t1 = has_time[2]
      }
      out = crop(out, extent = ee)
    }
  }
  return(out)
}

  