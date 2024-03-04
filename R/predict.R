#' Model prediction
#'
#' Apply a trained model on all pixels of a data cube.
#'
#' @param object a data cube proxy object (class cube)
#' @param model model used for prediction (e.g. from \code{caret} or \code{tidymodels})
#' @param ... further arguments passed to the model-specific predict method
#' @param output_names optional character vector for output variable(s)
#' @param keep_bands logical; keep bands of input data cube, defaults to FALSE, i.e. original bands will be dropped
#' @details 
#' The model-specific predict method will be automatically chosen based on the class of the provided model. It aims at supporting 
#' models from the packages \code{tidymodels}, \code{caret}, and simple models as from \code{lm} or \code{glm}. 
#' 
#' 
#' For multiple output variables or output in form of lists or data.frames, \code{output_names} must be provided and match 
#' names of the columns / items of the result object returned from the underlying predict method. For example, 
#' predictions using \code{tidymodels} return a tibble (data.frame) with columns like \code{.pred_class} (classification case).
#' This must be explicitly provided as \code{output_names}. Similarly, \code{predict.lm} and the like return lists
#' if the standard error is requested by the user and \code{output_names} hence should be set to \code{c("fit","se.fit")}.
#' 
#' For more complex cases or when predict expects something else than a \code{data.frame}, this function may not work at all. 
#' @note This function returns a proxy object, i.e., it will not immediately start any computations. 
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
predict.cube <- function(object, model, ..., output_names=c("pred"), keep_bands=FALSE) {
  
  arguments <- list(...)
  
  srcfile2 =  tempfile(pattern=".stream_", fileext = ".R")
  srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix
  
  # support custom library paths
  cat(paste0(".libPaths(",  paste(deparse(.libPaths()),collapse=""), ")\n"), file = srcfile2, append = TRUE) 
  cat("require(gdalcubes)", "\n", file = srcfile2, append = TRUE)
  load_pkgs = .packages()

  cat(paste0("require(", load_pkgs,")",collapse  = "\n"), "\n", file = srcfile2, append = TRUE) 
  
  # create new environment for model and args
  env = new.env(parent = .GlobalEnv)
  env$model = model
  env$args = arguments
  env$output_names = output_names
  
  load_env = env
  
  # if (sum(sapply(ls(envir = load_env), FUN = function(x) {object.size(get(x, envir = load_env))})) > 100*1024^2) {
  #   warning("The current environment seems to be rather large (> 100 Mb), this may result in reduced performance.")
  # }
  envfile = tempfile(pattern="renv_", fileext = ".rda")
  envfile = gsub("\\\\", "/", envfile) # Windows fix
  save(list = ls(envir = load_env),file = envfile, envir = load_env)
  cat(paste0("load(\"", envfile, "\")"), "\n", file = srcfile2, append = TRUE)
  script = system.file("scripts/f_predict.R",package = "gdalcubes")
  cat(paste("eval(parse(\"", script, "\"))", sep=""), "\n", file = srcfile2, append = TRUE)
  cmd <- paste(file.path(R.home("bin"),"Rscript"), " --vanilla ", srcfile2, sep="")
  
  # for tidymodels, support .pred columns without warning
  output_names = sapply(output_names, function(s) {
    if (!grepl("^[[:alnum:]].*", s)) {
      return(paste0("X",s))
    }
    return(s)
  })
  
  x = gc_create_stream_apply_pixel_cube(object, cmd, length(output_names), output_names, keep_bands)
  class(x) <- c("apply_pixel_cube", "cube", "xptr")
  return(x)
  
}