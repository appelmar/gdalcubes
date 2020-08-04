#' Convert complete image collections to cloud-optimized GeoTIFFs
#'
#' This function translates all images of a gdalcubes image collection to cloud-optimized GeoTIFF files.
#' The output contains converted imagery as an additional copy (original files are not deleted) and
#' a new image collection file.
#' 
#' @param collection path to an existing image collection file
#' @param target_dir directory where the output will be stored, will be created if necessary
#' @param overwrite logical; if TRUE existing files will be overwritten
#' @param creation_options further settings of the GDAL COG driver; see \url{https://gdal.org/drivers/raster/cog.html}
#' @return path to the new image collection file for use as argument to \code{\link{image_collection}} 
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
#' L8.col
#' \donttest{
#' if ("COG" %in% gdalcubes_gdalformats()) {
#'   L8.cog.col = translate_cog(L8.col)
#'   L8.cog.col
#' }
#' }
#' @note This function requires the GDAL COG driver, which was added in GDAL version 3.1.
#' @export
translate_cog <- function(collection, target_dir = tempfile(pattern = "image_collection_"), overwrite = TRUE, creation_options = c("BLOCKSIZE=256","COMPRESS=DEFLTE","LEVEL=1","RESAMPLING=CUBIC")) {
  stopifnot(is.image_collection(collection))
  out = libgdalcubes_translate_cog(collection, target_dir, .pkgenv$threads, overwrite, creation_options)
  return(out)
}




#' Convert complete image collections to cloud-optimized GeoTIFFs
#'
#' This function translates all images of a gdalcubes image collection to GeoTIFF files.
#' The output contains converted imagery as an additional copy (original files are not deleted) and
#' a new image collection file.
#' 
#' @param collection path to an existing image collection file
#' @param target_dir directory where the output will be stored, will be created if necessary
#' @param overwrite logical; if TRUE existing files will be overwritten
#' @param creation_options further settings of the GDAL GTiff driver; see \url{https://gdal.org/drivers/raster/gtiff.html}
#' @return path to the new image collection file for use as argument to \code{\link{image_collection}} 
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
#' L8.col
#' 
#' \donttest{
#' L8.tif.col = translate_gtiff(L8.col)
#' L8.tif.col
#' }
#' 
#' @details 
#' The functions \code{\link{translate_gtiff}} and \code{\link{translate_cog}} have the same purpose to convert imagery to optimized GeoTIFF files. The latter uses the recent COG GDAL driver,
#' whereas the former uses the normal GTiff driver. Depending on additional creation options and the input images, files creted with \code{\link{translate_gtiff}} may or may not contain overview images.
#' @export
translate_gtiff <- function(collection, target_dir = tempfile(pattern = "image_collection_"), overwrite = TRUE, creation_options = c("TILED=YES","COMPRESS=DEFLATE","ZLEVEL=1","COPY_SRC_OVERVIEWS=TRUE")) {
  stopifnot(is.image_collection(collection))
  out = libgdalcubes_translate_gtiff(collection, target_dir, .pkgenv$threads, overwrite, creation_options)
  return(out)
}

