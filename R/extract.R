#' Extract values from a data cube by spatial or spatiotemporal features
#'
#' Extract pixel values of a data cube from a set of spatial or spatiotemporal features. 
#' Applications include the extraction of full time 
#' series at irregular points, extraction from spatiotemporal points, extraction of
#' pixel values in polygons, and computing summary statistics over polygons.
#' 
#' @param cube source data cube to extract values from
#' @param sf object of class \code{sf}, see \link[sf:st_as_sf]{sf package}
#' @param datetime Date, POSIXt, or character vector containing per feature time information; length must be identical to the number of features in \code{sf}
#' @param time_column name of the column in \code{sf} containing per feature time information
#' @param FUN optional function to compute per feature summary statistics
#' @param ... additional arguments passed to \code{FUN}
#' @param reduce_time logical; if TRUE, time is ignored when \code{FUN} is applied
#' @return A data.frame with columns FID, time, and data cube bands / variables, see Details 
#' @details 
#' The geometry in \code{sf} can be of any simple feature type supported by GDAL, including 
#' POINTS, LINES, POLYGONS, MULTI*, and more. If no time information is provided
#' in one of the arguments \code{datetime} or \code{time_column}, the full time series
#' of pixels with regard to the features are returned. 
#' 
#' Notice that feature identifiers in the \code{FID} column typically correspond to the row names / numbers 
#' of the provided sf object. This can be used to combine the output with the original geometries, e.g., using \code{\link[base:merge]{merge()}}.
#'  
#' Pixels with missing values are automatically dropped from the result. It is hence not
#' guaranteed that the result will contain rows for all input features.
#'
#' Features are automatically reprojected if the coordinate reference system differs from the data cube.
#' 
#' Extracted values can be aggregated by features by providing a summary function. 
#' If \code{reduce_time} is FALSE (the default), the values are grouped 
#' by feature and time, i.e., the result will contain unique combinations of FID and time.
#' To ignore time and produce a single value per feature, \code{reduce_time} can be set to TRUE.
#' 
#' @examples
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE)
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(srs="EPSG:32618", dy=1000, dx=1000, dt="P1M",
#'               aggregation = "median", resampling = "bilinear",
#'               extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931,
#'                           t0="2018-01-01", t1="2018-04-30"))
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube = select_bands(L8.cube, c("B04", "B05"))
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI")
#' L8.ndvi
#' 
#' if (gdalcubes_gdal_has_geos()) {
#'   if (requireNamespace("sf", quietly = TRUE)) {
#'   
#'     # create 50 random point locations
#'     x = runif(50, v$space$left, v$space$right)
#'     y = runif(50, v$space$bottom, v$space$top)
#'     t = sample(seq(as.Date("2018-01-01"),as.Date("2018-04-30"), by = 1),50, replace = TRUE)
#'     df = sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = v$space$srs)
#' 
#'     # 1. spatiotemporal points
#'     extract_geom(L8.ndvi, df, datetime = t)
#' 
#'     \donttest{
#'     # 2. time series at spatial points
#'     extract_geom(L8.ndvi, df)
#'   
#'     # 3. summary statistics over polygons
#'     x = sf::st_read(system.file("nycd.gpkg", package = "gdalcubes"))
#'     zstats = extract_geom(L8.ndvi,x, FUN=median, reduce_time = TRUE)
#'     zstats
#'     # combine with original sf object
#'     x$FID = rownames(x)
#'     x = merge(x, zstats, by = "FID")
#'     x
#'     plot(x["NDVI"])
#'     }
#'   }
#' }
#' @export
extract_geom = function(cube, sf, datetime = NULL, time_column = NULL, FUN = NULL, ..., reduce_time = FALSE) {
  
  stopifnot(is.cube(cube))
  if (!is.null(datetime) && !is.null(time_column)) {
    warning("expected only one of datetime or time_column; time_column will be ignored")
    time_column = NULL
  }
  if (!is.null(FUN)) {
    if (!is.function(FUN)) {
      stop("FUN is not a function")
    }
  }
  
  stopifnot("sf" %in% class(sf)) 
  if (!requireNamespace("sf", quietly = TRUE))
    stop("sf package not found, please install first") 
  
  
  if (is.null(time_column) && !is.null(datetime)) {
    if (length(datetime) != nrow(sf)) {
      stop("Length of datetime does not match the number of features")
    }
    time_column = basename(tempfile(pattern="time_")) # random name with prefix "time_"
    if (!is.character(datetime)) {
      datetime = as.character(datetime)
    }
    # add column
    sf[time_column] = datetime
  }
  else if (is.null(time_column) && is.null(datetime)) {
    time_column = ""
  }
  
  ogrfile = tempfile(fileext = ".gpkg")
  sf::st_write(sf, ogrfile, quiet = TRUE)
  
  out = gc_extract(cube, ogrfile, time_column)
  if (!is.null(FUN)) {
    if (!reduce_time) {
      out = stats::aggregate(out[,-(1:2)], by = list(out$FID, out$time), FUN, ...)
      names(out) = c("FID", "time", names(cube))
    }
    else {
      out = stats::aggregate(out[,-(1:2)], by = list(out$FID), FUN, ...)
      names(out) = c("FID", names(cube))
    }
    
  }

  # remove temporary ogrfile if still exists
  if (file.exists(ogrfile)) {
    file.remove(ogrfile)
  }
  
  return(out)
}

  
  
  
