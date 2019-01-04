


#' 
#' Create a spatiotemporal data cube view
#' 
#' Data cube views define the shape of a cube, i.e., the spatiotemporal extent, resolution, and projection how to look at the data.
#' They are used to access image collections as on-demand data cubes. The data cube will filter images based on the view's
#' extent, read image data at the defined resolution, and warp to the output projection automatically. All parameters
#' are optional. By default gdalcubes tries to derive a view that covers the whole image collection at course resolution, which
#' may or may not make sense in practice. It is recommended to at least spcify the spatial and temporal extent, the resolution (nx, ny, dt) 
#' and the output projection.
#' 
#' @note All arguments are optional
#' 
#' @param cube Return the view of the specified data cube, other arguments will be ignored
#' @param view if provided, update the view object instead of creating a new data cube view where fields that are already set will be overwritten
#' @param proj output projection as string, can be proj4 definition, WKT, or in the form "EPSG:XXXX"
#' @param nx number of pixels in x-direction (longitude / easting)
#' @param ny number of pixels in y-direction (latitude / northing)
#' @param dx size of pixels in x-direction (longitude / easting)
#' @param dy size of pixels  in y-direction (latitude / northing)
#' @param l left boundary of the spatial extent (expressed in the coordinate system defined by proj or WGS84 by default)
#' @param r right boundary of the spatial extent (expressed in the coordinate system defined by proj or WGS84 by default)
#' @param t top boundary of the spatial extent (expressed in the coordinate system defined by proj or WGS84 by default)
#' @param b bottom boundary of the spatial extent (expressed in the coordinate system defined by proj or WGS84 by default)
#' @param dt size of pixels in time-direction, expressed as ISO8601 period string (only 1 number and unit is allowed) such as "P16D"
#' @param t0 start date/time of the temporal extent (expressed as ISO8601 datetime string)
#' @param nt number of pixels in t-direction
#' @param t1 end date/time of the temporal extent (expressed as ISO8601 datetime string)
#' @param aggregation aggregation method as string, defining how to deal with pixels with data from multiple images, can be "min", "max", "mean", "median", "first"
#' @param resampling resampling method used in gdalwarp when images are read, can be "near", "bilinear", "bicubic" or others as supported by gdalwarp (see \url{https://www.gdal.org/gdalwarp.html})
#' @examples 
#' gcbs_view(proj="EPSG:4326", l = -20, r = 20, t = 60, b=40, 
#'           t0="2018-01-01", t1="2018-09-30", dt="P1M", nx=1000, 
#'           ny=500, aggregation = "mean", resampling="bilinear")
#' @return A list with view properties
#' @export
gcbs_view <- function(cube, view, proj, nx, ny, dx, dy, l, r, t, b, dt, nt, t0, t1, aggregation,resampling)  {
  if (!missing(cube)) {
    if (length(as.list(match.call())) > 2) {
      warning("Provided arguments except cube will be ignored")
    }
    stopifnot(is.gcbs_cube(cube))
    x = libgdalcubes_get_cube_view(cube)
    class(x) <- c("gcbs_view", class(x))
    return(x)
  }
  xx = NULL
  if (missing(view)) {
    xx = list(space =
                list(left = NULL,
                     right = NULL,
                     top = NULL,
                     bottom = NULL,
                     nx = NULL,
                     dx = NULL,
                     ny = NULL,
                     dy = NULL,
                     proj=NULL),
              time = list(
                t0 = NULL,
                t1 = NULL,
                dt = NULL,
                nt = NULL
              ),
              aggregation = NULL,
              resampling = NULL
              )
  }
  else {
    stopifnot(is.gcbs_view(view))
    xx = view
  }

  if (!(missing(l) && missing(r) && missing(t) && missing(b)) && missing(proj)) {
    warning("Spatial extent without coordinate system given, assuming WGS84")
    proj = "EPSG:4326"
  }
  
  if (!missing(proj)) xx$space$proj = proj
  if (!missing(nx)) xx$space$nx = as.integer(nx)
  if (!missing(ny)) xx$space$ny = as.integer(ny)
  if (!missing(dx)) xx$space$dx = as.double(dx)
  if (!missing(dy)) xx$space$dy = as.double(dy)
  if (!missing(l)) xx$space$left = as.double(l)
  if (!missing(r)) xx$space$right = as.double(r)
  if (!missing(b)) xx$space$bottom = as.double(b)
  if (!missing(t)) xx$space$top = as.double(t)
  if (!missing(t0)) {
    if (is.character(t0)) {
      xx$time$t0 = as.character(t0)
    }
    else if ("POSIXt" %in% class(t0)) {
      xx$time$t0 = format(t0, "%Y-%m-%dT%H:%M:%S")
    }
    else if ("Date" %in% class(t0)) {
      xx$time$t0 = format(t0, "%Y-%m-%d")
    }
    else {
      warning("Invalid type for t0, expected ISO 8601 character, POSIXt, or Date object, value will be ignored")
    }
  }
  if (!missing(t1)) {
    if (is.character(t1)) {
      xx$time$t1 = as.character(t1)
    }
    else if ("POSIXt" %in% class(t1)) {
      xx$time$t1 = format(t1, "%Y-%m-%dT%H:%M:%S")
    }
    else if ("Date" %in% class(t1)) {
      xx$time$t1 = format(t1, "%Y-%m-%d")
    }
    else {
      warning("Invalid type for t1, expected ISO 8601 character, POSIXt, or Date object, value will be ignored")
    }
  }
  if (!missing(dt)) {
    dt = toupper(as.character(dt))
    if (substr(dt[1],1,1) != "P") {
      dt[1] = paste("P", dt[1], sep="")
    }
    xx$time$dt = dt[1]
  }
  if (!missing(nt)) xx$time$nt = as.integer(nt)
  
  if (!missing(aggregation)) xx$aggregation = aggregation
  if (!missing(resampling)) xx$resampling = resampling
  
  class(xx) <- c("gcbs_view", class(xx))
  return(xx)
}



is.gcbs_view <- function(obj) {
 return("gcbs_view" %in% class(obj))
}

