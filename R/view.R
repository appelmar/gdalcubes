
#' Create or update a spatiotemporal data cube view
#'
#' Data cube views define the shape of a cube, i.e., the spatiotemporal extent, resolution, and spatial reference system (srs).
#' They are used to access image collections as on-demand data cubes. The data cube will filter images based on the view's
#' extent, read image data at the defined resolution, and warp / reproject images to the target srs automatically. 
#' 
#' @param view if provided, update this cube_view object instead of creating a new data cube view where fields that are already set will be overwritten
#' @param extent spatioptemporal extent as a list e.g. from \code{\link{extent}} or an image collection object, see Details
#' @param srs target spatial reference system as a string; can be a proj4 definition, WKT, or in the form "EPSG:XXXX"
#' @param nx number of pixels in x-direction (longitude / easting)
#' @param ny number of pixels in y-direction (latitude / northing)
#' @param dx size of pixels in x-direction (longitude / easting)
#' @param dy size of pixels in y-direction (latitude / northing)
#' @param dt size of pixels in time-direction, expressed as ISO8601 period string (only 1 number and unit is allowed) such as "P16D"
#' @param nt number of pixels in t-direction
#' @param keep.asp if TRUE, derive ny or dy automatically from nx or dx (or vice versa) based on the aspect ratio of the spatial extent
#' @param aggregation aggregation method as string, defining how to deal with pixels containing data from multiple images, can be "min", "max", "mean", "median", or "first"
#' @param resampling resampling method used in gdalwarp when images are read, can be "near", "bilinear", "bicubic" or others as supported by gdalwarp (see \url{https://gdal.org/programs/gdalwarp.html})
#' @details 
#' The \code{extent} argument expects a simple list with elementes \code{left}, \code{right}, \code{bottom}, \code{top}, \code{t0} (start date/time), \code{t1} (end date/time) or an image collection object.
#' In the latter case, the \code{\link{extent}} function is automatically called on the image collection object to get the full spatiotemporal extent of the collection. In the former case, datetimes 
#' are expressed as ISO8601 datetime strings.
#' 
#' The function can be used in two different ways. First, it can create data cube views from scratch by defining the extent, the spatial reference system, and for each dimension either the cell size (dx, dy, dt) or the total number of cells (nx, ny, nt).
#' Second, the function can update an existing data cube view by overwriting specific fields. In this case, the extent or some elements of the extent may be missing. 
#'
#' In some cases, the extent of the view is automatically extended if the provided resolution would end within a pixel. For example,
#' if the spatial extent covers an area of 1km x 1km and dx = dy = 300m, the extent would be enlarged to 1.2 km x 1.2km. The alignment will be reported to the user in 
#' a diagnostic message.  
#'
#'
#' @examples
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  L8.col = create_image_collection(L8_files, "L8_L1TP")
#'  
#'  # 1. Create a new data cube view specification
#'  v = cube_view(extent=extent(L8.col,"EPSG:4326"), srs="EPSG:4326", dt="P1M", 
#'            nx=1000, ny=500, aggregation = "mean", resampling="bilinear")
#'  v
#'
#'  # 2. overwrite parts of an existing data cube view
#'  vnew = cube_view(v, dt="P1M")
#' @return A list with data cube view properties
#' @export
cube_view <- function(view, extent, srs, nx, ny, nt, dx, dy, dt, aggregation, resampling, keep.asp=TRUE) {
  
  
  # general parameter checks
  if (!missing(view)) 
    stopifnot(is.cube_view(view))
  if (!missing(extent)) stopifnot(is.list(extent) || is.image_collection(extent))
  if (!missing(nx)) stopifnot(length(nx) == 1, nx %% 1 == 0)
  if (!missing(ny)) stopifnot(length(ny) == 1, ny %% 1 == 0)
  if (!missing(nt)) stopifnot(length(nt) == 1, nt %% 1 == 0)
  if (!missing(dx)) stopifnot(length(dx) == 1)
  if (!missing(dy)) stopifnot(length(dy) == 1)
  if (!missing(dt)) stopifnot(length(dt) == 1)
  if (!missing(aggregation)) stopifnot(is.character(aggregation),  length(aggregation) == 1)
  if (!missing(resampling)) stopifnot(is.character(resampling),  length(resampling) == 1)
  if (!missing(srs)) stopifnot(is.character(srs),  length(srs) == 1)

  xx = list(space =
              list(left = NULL,
                   right = NULL,
                   top = NULL,
                   bottom = NULL,
                   nx = NULL,
                   dx = NULL,
                   ny = NULL,
                   dy = NULL,
                   srs=NULL),
            time = list(
              t0 = NULL,
              t1 = NULL,
              dt = NULL,
              nt = NULL
            ),
            aggregation = NULL,
            resampling = NULL)
  
  
  
  if (!missing(view)) {
    if (missing(srs)) {
      srs = view$space$srs
    }
    xx$space$srs = srs
    
    # update existing view
    extent_old=list(left   = view$space$left,
                    right  = view$space$right,
                    top    = view$space$top,
                    bottom = view$space$bottom,
                    t0     = view$time$t0,
                    t1     = view$time$t1)
    
     # update extent first
     if (missing(extent)) {
       extent = extent_old
     }
    else {
      if (is.image_collection(extent)) {
        extent = extent(extent, srs) # wow!
      }
      else {
        if (is.null(extent$left)) extent$left = extent_old$left
        if (is.null(extent$right)) extent$right = extent_old$right
        if (is.null(extent$top)) extent$top = extent_old$top
        if (is.null(extent$bottom)) extent$bottom = extent_old$bottom
        if (is.null(extent$t0)) extent$t0 = extent_old$t0
        if (is.null(extent$t1)) extent$t1 = extent_old$t1
      }
    }
    
    xx$space$left = extent$left
    xx$space$right = extent$right
    xx$space$top = extent$top
    xx$space$bottom = extent$bottom
    xx$time$t0 = extent$t0
    xx$time$t1 = extent$t1
    # now we have a well-defined extent (though we still should check for meaningful values here)
    
    
    
    if (!missing(nx) && !missing(dx)) {
      warning("conflicting arguments nx and dx, ignoring dx")
    }
    if (!missing(ny) && !missing(dy)) {
      warning("conflicting arguments ny and dy, ignoring dy")
    }
    if (!missing(dt) && !missing(nt)) {
      warning("conflicting arguments nt and dt, ignoring nt")
    }

    xres_defined = FALSE
    yres_defined = FALSE
    
    if (!missing(nx) || !missing(dx)) {
      if (missing(nx)) {
        xx$space$dx = dx
      }
      else {
        xx$space$nx = nx
      }
      xres_defined = TRUE
    }
    
    if (!missing(ny) || !missing(dy)) {
      if (missing(ny)) {
        xx$space$dy = dy
      }
      else {
        xx$space$ny = ny
      }
      yres_defined = TRUE
    }
    
    # try to derive x resolution from y by keeping the aspcet ratio of the extent
    if (!xres_defined) {
      if (keep.asp) {
        if (!missing(ny)) {
          xx$space$nx = round(ny * (extent$right - extent$left)/(extent$top - extent$bottom) )
          xres_defined = TRUE  
        }
        else if (!missing(dy)) {
          xx$space$dx = dy
          xres_defined = TRUE
        }
      }
    }
    
    if (!yres_defined) {
      if (keep.asp) {
        if (!missing(nx)) {
          xx$space$ny = round(nx * (extent$top - extent$bottom)/(extent$right - extent$left))
          yres_defined = TRUE  
        }
        else if (!missing(dx)) {
          xx$space$dy = dx
          yres_defined = TRUE
        }
      }
    }
    
    # use from input view
    if (!xres_defined) {
       if (!is.null(view$space$nx)) {
         xx$space$nx = view$space$nx
         xres_defined = TRUE 
       } 
      else if (!is.null(view$space$dx)) {
        xx$space$dx = view$space$dx
        xres_defined = TRUE 
      } 
    }
    if (!yres_defined) {
      if (!is.null(view$space$ny)) {
        xx$space$ny = view$space$ny
        yres_defined = TRUE 
      } 
      else if (!is.null(view$space$dy)) {
        xx$space$dy = view$space$dy
        yres_defined = TRUE 
      } 
    }
    
    
    if (! xres_defined) {
      stop("definition of x dimension is incomplete")
    }
    if (! yres_defined) {
      stop("definition of y dimension is incomplete")
    }
    
    
    
    if (!missing(nt) || !missing(dt)) {
      if (missing(dt)) {
        xx$time$nt = nt
      }
      else {
        xx$time$dt = dt
      }
    }
    else {
      if (!is.null(view$time$dt)) {
        xx$time$dt = view$time$dt
      }
      else if (!is.null(view$time$nt)) {
        xx$time$nt = view$time$nt
      }
    }
    
    xx$aggregation = ifelse(missing(aggregation), view$aggregation, aggregation)
    xx$resampling = ifelse(missing(resampling), view$resampling, resampling)
  

  }
  else {
    # create a new data cube view from scratch
    
    if (missing(srs)) {
      warning("argument srs is missing, assuming WGS84")
      srs = "EPSG:4326"
    }
    xx$space$srs = srs
    
    if (missing(extent)) {
      stop("argument extent is required")
    }
    
    if (is.image_collection(extent)) {
      extent = extent(extent, srs) # wow!
    }
    else {
      if (is.null(extent$left)) {
        stop("argument extent does not contain left boundary")
      }
      if (is.null(extent$right)) {
        stop("argument extent does not contain right boundary")
      }
      if (is.null(extent$bottom)) {
        stop("argument extent does not contain bottom boundary")
      }
      if (is.null(extent$top)) {
        stop("argument extent does not contain top boundary")
      }
      if (is.null(extent$t0)) {
        stop("argument extent does not contain t0 (start date/time) boundary")
      }
      if (is.null(extent$t1)) {
        stop("argument extent does not contain t1 (end date/time) boundary")
      }
    }
    xx$space$left = extent$left
    xx$space$right = extent$right
    xx$space$top = extent$top
    xx$space$bottom = extent$bottom
    xx$time$t0 = extent$t0
    xx$time$t1 = extent$t1
    
    # now we have a well-defined extent (though we still should check for meaningful values here)
    
  
    
    if (!missing(nx) && !missing(dx)) {
      warning("conflicting arguments nx and dx, ignoring dx")
    }
    if (!missing(ny) && !missing(dy)) {
      warning("conflicting arguments ny and dy, ignoring dy")
    }
    if (!missing(dt) && !missing(nt)) {
      warning("conflicting arguments nt and dt, ignoring nt")
    }
    
    xres_defined = FALSE
    yres_defined = FALSE
    
    if (!missing(nx) || !missing(dx)) {
      if (missing(nx)) {
        xx$space$dx = dx
      }
      else {
        xx$space$nx = nx
      }
      xres_defined = TRUE
    }
    
    if (!missing(ny) || !missing(dy)) {
      if (missing(ny)) {
        xx$space$dy = dy
      }
      else {
        xx$space$ny = ny
      }
      yres_defined = TRUE
    }
    
    # try to derive x resolution from y by keeping the aspcet ratio of the extent
    if (!xres_defined) {
      if (keep.asp) {
        if (!missing(ny)) {
          xx$space$nx = round(ny * (extent$right - extent$left)/(extent$top - extent$bottom) )
          xres_defined = TRUE  
        }
        else if (!missing(dy)) {
          xx$space$dx = dy
          xres_defined = TRUE
        }
      }
    }
    
    if (!yres_defined) {
      if (keep.asp) {
        if (!missing(nx)) {
          xx$space$ny = round(nx * (extent$top - extent$bottom)/(extent$right - extent$left))
          yres_defined = TRUE  
        }
        else if (!missing(dx)) {
          xx$space$dy = dx
          yres_defined = TRUE
        }
      }
    }
    
    if (! xres_defined) {
      stop("definition of x dimension is incomplete, one of nx and dx is required")
    }
    if (! yres_defined) {
      stop("definition of y dimension is incomplete, one of ny and dy is required")
    }
    
    
    if (!missing(nt) || !missing(dt)) {
      if (missing(dt)) {
        xx$time$nt = nt
      }
      else {
        xx$time$dt = dt
      }
    }
    else {
      stop("definition of t dimension is incomplete, one of nt and dt is required")
    }
    
    xx$aggregation = ifelse(missing(aggregation), "first", aggregation)
    xx$resampling = ifelse(missing(resampling), "nearest", resampling)
    
  }
  
  
  # 
  xx = libgdalcubes_create_view(xx)
  
  class(xx) <- c("cube_view", class(xx))
  return(xx)

}


is.cube_view <- function(obj) {
  return("cube_view" %in% class(obj))
}

#' Print data cube view information
#' 
#' Prints information about a data cube view, including its dimensions, spatial reference, aggregation method, and resampling method.
#' 
#' @param x Object of class "cube_view"
#' @param ... Further arguments passed to the generic print function
#' @examples 
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' print(v)
#' @export
print.cube_view <- function(x, ...) {
  stopifnot(is.cube_view(x))
  cat("A data cube view object\n\n")
  cat("Dimensions:\n")
  dims = data.frame(low   = c(x$space$left, x$space$bottom, x$time$t0),
                    high  = c(x$space$right, x$space$top, x$time$t1),
                    count = c(x$space$nx, x$space$ny, x$time$nt),
                    pixel_size  = c(x$space$dx, x$space$dy, x$time$dt))
  rownames(dims) <- c("x", "y", "t")
  dims = dims[3:1, ] # reverse rows
  print(dims)
  cat(paste("\nSRS: \"", x$space$srs, "\"\n", sep=""))
  if (!is.null(x$aggregation))
    cat(paste("Temporal aggregation method: \"", x$aggregation, "\"\n", sep=""))
  if(!is.null(x$resampling))
    cat(paste("Spatial resampling method: \"", x$resampling, "\"\n", sep=""))
  cat("\n")
}
