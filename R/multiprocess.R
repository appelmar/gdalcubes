
#' Internal function calld in a gdalcubes worker process
#' 
#' 
#' @noRd
.exec_worker <- function() {
  if(!Sys.getenv("GDALCUBES_WORKER") == "1") {
    stop("Function can only be called within a worker process")
  }
  
  json = Sys.getenv("GDALCUBES_WORKER_JSON")
  stopifnot(file.exists(json))
  
  worker_id = as.integer(Sys.getenv("GDALCUBES_WORKER_ID"))
  stopifnot(is.finite(worker_id) && worker_id >= 0)
  
  worker_n = as.integer(Sys.getenv("GDALCUBES_WORKER_N"))
  stopifnot(is.finite(worker_n) && worker_n > 0)
  
  work_dir = Sys.getenv("GDALCUBES_WORKER_DIR")
  stopifnot(dir.exists(work_dir))
  
  gdalcubes_options(log_file = file.path(work_dir, paste0(sprintf("worker_%03d", worker_id), ".log")))
  
  gc_exec_worker(json, worker_id, worker_n, work_dir)
  
}







