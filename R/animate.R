

#' Animate a data cube as an image time series
#'
#' @param x a data cube proxy object (class cube)
#' @param ... parameters passed to plot.cube
#' @param fps frames per second of the animation
#' @param loop how many iterations, 0 = infinite
#' @param width width (in pixels) of the animation
#' @param height height (in pixels) of the animation
#' @param save_as character path where the animation shall be stored as a gif file
#' @param plot logical; plot the animation (default is TRUE) 
#' @details 
#' Animations can be created for single band data cubes or RGB plots of multi-band data cubes (by providing the argument rgb) only.
#' @seealso \code{\link{plot.cube}}
#' @examples 
#' \donttest{
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"))
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P16D")
#' 
#' animate(select_bands(raster_cube(L8.col, v), c("B02", "B03", "B04")), rgb=3:1,
#'         zlim=c(0,20000), fps=1, loop=1)
#' 
#' animate(select_bands(raster_cube(L8.col, v), c("B05")), col=terrain.colors, key.pos=1)
#' }
#' @export
animate  <- 
  function(x,
           ...,
           fps = 1,
           loop = 0,
           width = dev.size(units = "px")[1],
           height = dev.size(units = "px")[2],
           save_as = NULL,
           plot = TRUE
           ) {
    
    if (!requireNamespace("magick", quietly = TRUE))
      stop("magick package not found, please install first") 
    
    if(is.null(save_as) && !plot) {
      stop("nothing to do, please set either plot = TRUE or save_as to a filename")
    }
    
    stopifnot(is.cube(x))
    size = c(nbands(x), size(x))
    
    if(size[2] == 1) {
      stop("nothing to animate; cube has only one time slice")
    }
    
    additional_args = list(...)
    if (is.null(additional_args$rgb) && size[1] > 1) {
      stop("animate works only for RGB plots, or single band data cubes")
    }
      
  
    fname_start = tempfile()
    png(filename = paste(fname_start, "_%04d.png", sep=""), width=width, height=height)

    tryCatch({
      for (i in 1:size[2]) {
        args = additional_args
        args$t = i
        args$x = x
        do.call("plot.cube", args=args)
      }}, finally = {dev.off()})
    
    imgs = list.files(dirname(fname_start), pattern = paste(basename(fname_start), ".*\\.png" , sep=""), full.names = TRUE)
    frames = magick::image_read(imgs[1])
    for (i in 2:length(imgs)) {
      frames = c(frames,  magick::image_read(imgs[i]))
    }
    animation = magick::image_animate(frames, fps=fps, loop = loop)
    if (!is.null(save_as)) {
      magick::image_write(animation, save_as)
    }
    if (plot)
      print(animation)
    return(animation)
  }
