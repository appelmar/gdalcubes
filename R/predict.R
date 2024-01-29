#' Predict values based on a trained model to all pixels of a data cube.
#'
#' @param object a data cube proxy object (class cube)
#' @param model model used for prediction (e.g. from \code{caret})
#' @param ... further arguments passed to predict
#' @details 
#' The predict method used will be automatically chosen based on the model class. It should support models from the \code{caret} 
#' package and simple models as from \code{lm} or \code{glm} . However, it  currently works only for a single output variable and 
#' requires a \code{data.frame} input. 
#' @examples 
#' \donttest{
#' # create image collection from example Landsat data only
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE)
#' }
#' 
#' v = cube_view(extent=list(left=388941.2, right=766552.4,
#'                           bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P3M")
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' 
#' x = sf::st_read(system.file("ny_samples.gpkg", package = "gdalcubes"))
#' 
#' raster_cube(L8.col, v) |>
#'   select_bands(c("B02","B03","B04","B05")) |>
#'   extract_geom(x) -> train
#' 
#' x$FID = rownames(x)  
#' train = merge(train, x, by = "FID")
#' train$iswater = as.factor(train$class == "water")
#' 
#' log_model <- glm(iswater ~ B02 + B03 + B04 + B05, data = train, family = "binomial")
#' 
#' raster_cube(L8.col, v) |>
#'   select_bands(c("B02","B03","B04","B05")) |>
#'   predict(model=log_model, type="response") |>
#'   plot(key.pos=1)
#' }
#' @export
predict.cube <- function(object, model, ...) {
  
  arguments <- list(...)
  paste(arguments)
  
  srcfile2 =  tempfile(pattern=".stream_", fileext = ".R")
  srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix
  
  # Create new script, starting with conditional initialization in order to
  # avoid reloading packages / environment for every chunk processed in the current process
  cat("if(!exists(\"IS_INITIALIZED\") || IS_INITIALIZED==FALSE) {", "\n", file = srcfile2, append = FALSE) # if
  
  # support custom library paths
  cat(paste0(".libPaths(",  paste(deparse(.libPaths()),collapse=""), ")\n"), file = srcfile2, append = TRUE) 
  cat("require(gdalcubes)", "\n", file = srcfile2, append = TRUE)
  load_pkgs = .packages()

  cat(paste0("require(", load_pkgs,")",collapse  = "\n"), "\n", file = srcfile2, append = TRUE) 
  
  # create new environment for model and args
  env = new.env(parent = .GlobalEnv)
  env$model = model
  env$args = arguments
  
  load_env = env
  
  # if (sum(sapply(ls(envir = load_env), FUN = function(x) {object.size(get(x, envir = load_env))})) > 100*1024^2) {
  #   warning("The current environment seems to be rather large (> 100 Mb), this may result in reduced performance.")
  # }
  envfile = tempfile(pattern="renv_", fileext = ".rda")
  save(list = ls(envir = load_env),file = envfile, envir = load_env)
  cat(paste0("load(\"", envfile, "\")"), "\n", file = srcfile2, append = TRUE)
  script = system.file("scripts/f_predict.R",package = "gdalcubes")
  cat("IS_INITIALIZED=TRUE","\n", "}", "\n", file = srcfile2, append = TRUE) # end if
  cat(paste("eval(parse(\"", script, "\"))", sep=""), "\n", file = srcfile2, append = TRUE)
  cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
  
  keep_bands = FALSE
  names = "prediction"
  nb = 1

  x = gc_create_stream_apply_pixel_cube(object, cmd, nb, names, keep_bands)
  class(x) <- c("apply_pixel_cube", "cube", "xptr")
  return(x)
  
}