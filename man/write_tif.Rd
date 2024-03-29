% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cube.R
\name{write_tif}
\alias{write_tif}
\title{Export a data cube as a collection of GeoTIFF files}
\usage{
write_tif(
  x,
  dir = tempfile(pattern = ""),
  prefix = basename(tempfile(pattern = "cube_")),
  overviews = FALSE,
  COG = FALSE,
  rsmpl_overview = "nearest",
  creation_options = NULL,
  write_json_descr = FALSE,
  pack = NULL
)
}
\arguments{
\item{x}{a data cube proxy object (class cube)}

\item{dir}{destination directory}

\item{prefix}{output file name}

\item{overviews}{logical; generate overview images}

\item{COG}{logical; create cloud-optimized GeoTIFF files (forces overviews=TRUE)}

\item{rsmpl_overview}{resampling method for overviews (image pyramid) generation (see \url{https://gdal.org/programs/gdaladdo.html} for available methods)}

\item{creation_options}{additional creation options for resulting GeoTIFF files, e.g. to define compression (see \url{https://gdal.org/drivers/raster/gtiff.html#creation-options})}

\item{write_json_descr}{logical; write a JSON description of x as additional file}

\item{pack}{reduce output file size by packing values (see Details), defaults to no packing}
}
\value{
returns (invisibly) a vector of paths pointing to the created GeoTIFF files
}
\description{
This function will time slices of a data cube as GeoTIFF files
in a given directory.
}
\details{
If \code{write_json_descr} is TRUE, the function will write an additional file with name according to prefix (if not missing) or simply cube.json 
This file includes a serialized description of the input data cube, including all chained data cube operations.

Additional GDAL creation options for resulting GeoTIFF files must be passed as a named list of simple strings, where element names refer to the key. For example,
\code{creation_options = list("COMPRESS" = "DEFLATE", "ZLEVEL" = "5")} would enable deflate compression at level 5.

To reduce the size of created files, values can be packed by applying a scale factor and an offset value and using a smaller
integer data type for storage. The \code{pack} argument can be either NULL (the default), or a list with elements \code{type}, \code{scale}, \code{offset}, 
and \code{nodata}. \code{type} can be any of "uint8", "uint16" , "uint32", "int16", or "int32". \code{scale}, \code{offset}, and 
\code{nodata} must be numeric vectors with length one or length equal to the number of data cube bands (to use different values for different bands). 
The helper function  \code{\link{pack_minmax}} can be used to derive offset and scale values with maximum precision from minimum and maximum data values on
original scale.

If \code{overviews=TRUE}, the numbers of pixels are halved until the longer spatial dimensions counts less than 256 pixels.
Setting \code{COG=TRUE} automatically sets \code{overviews=TRUE}.
}
\examples{
# create image collection from example Landsat data only 
# if not already done in other examples
if (!file.exists(file.path(tempdir(), "L8.db"))) {
  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
                         ".TIF", recursive = TRUE, full.names = TRUE)
  create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
}

L8.col = image_collection(file.path(tempdir(), "L8.db"))
v = cube_view(extent=list(left=388941.2, right=766552.4, 
              bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
              srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
write_tif(select_bands(raster_cube(L8.col, v), c("B04", "B05")), dir=tempdir())
}
\seealso{
\code{\link{pack_minmax}}
}
