#' Predict values based on a trained model to all pixels of a data cube.
#'
#' @param object a data cube proxy object (class cube)
#' @param model model used for prediction (e.g. from \code{caret})
#' @param ... further arguments (not used)
#' @details 
#' The predict method used will be automatically chosen based on the model class. It should support models from the \code{caret} 
#' package and simple models as from\code{lm}. However, it  currently works only for a single output variable and 
#' requires a \code{data.frame} input. 
#' @examples 
#' # TODO
#' 
#' @export
predict.cube <- function(object, model, ...) {
  
  # 1. find predict method for model
  model_type = class(model)
  mcand = paste0("predict.", model_type)
  a = mget(mcand, ifnotfound = NA, inherit = TRUE)
  predict_fun = NA
  for (i in 1:length(mcand)) {
    predict_fun = a[[mcand[i]]]
    if (is.function(predict_fun)) {
      break
    }
  }
  
  # 2. write script
  srcfile2 =  tempfile(pattern=".stream_", fileext = ".R")
  srcfile2 = gsub("\\\\", "/", srcfile2) # Windows fix
  
  
  # TODO: ADD IF TO AVOID RELOADING EVERYTHING FOR EACH CHUNK
  
  # Create new script, starting with conditional initialization in order to
  # avoid reloading packages / environment for every chunk processed in the current process
  cat("if(!exists(\"IS_INITIALIZED\") || IS_INITIALIZED==FALSE) {", "\n", file = srcfile2, append = FALSE) # if
  
  # support custom library paths
  cat(paste0(".libPaths(",  paste(deparse(.libPaths()),collapse=""), ")\n"), file = srcfile2, append = TRUE) 
  cat("require(gdalcubes)", "\n", file = srcfile2, append = TRUE)
  load_pkgs = .packages()

  cat(paste0("require(", load_pkgs,")",collapse  = "\n"), "\n", file = srcfile2, append = TRUE) 
  load_env = .GlobalEnv
  
  if (sum(sapply(ls(envir = load_env), FUN = function(x) {object.size(get(x, envir = load_env))})) > 100*1024^2) {
    warning("The current environment seems to be rather large (> 100 Mb), this may result in reduced performance.")
  }
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