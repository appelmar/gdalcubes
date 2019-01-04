#' Start one or more gdalcubes_server background processes
#'
#' Starts one or more gdalcubes_server instance on this machine. Created processes are added to a global list of gdalcubes_server processes.
#' Simultaneously running processes must use different ports.
#' 
#' @param port port number(s) where gdalcubes server(s) will listen for inccoming requests
#' @param endpoint base path where the API sets up its endpoints
#' @param whitelist character vector with hosts that are allowed to connect to the server, if null, all incoming requests will be accepted
#' @param threads number of threads running to process parallel chunk read requests
#' @param n number of servers to start, if n > 1 and ports has length 1, port numbers will be increased automatically by one
#' @export
gcbs_start_server <- function(port=1111, endpoint = "/gdalcubes/api", whitelist=NULL, threads=1, n=1) {
  stopifnot(requireNamespace("processx", quietly = TRUE))

  if (n > 1 && (length(port) != 1 && length(port) != n)) {
    warning(paste("Expected either 1 or n ports, using", paste(port[1]:(port[1]+n-1), collapse=",")))
    port = port[1]:(port[1]+n-1)
  }
  else if (n > 1 && length(port) == 1) {
    port = (port[1]+n-1)
  }

  # TODO: check input arguments
  # TODO: add workdir argument
  # TODO: add whitelist
  env <- .pkgenv
  if (!exists("gdalcubes.server_processes", envir=env)) {
    env$gdalcubes.server_processes = list()
  }
  #env$gdalcubes.server_processes[[length(env$server_processes) + 1]] <- processx::process$new("gdalcubes_server", "-b", endpoint, "-p", port, "-t", threads)


  for (i in 1:n) {
    env$gdalcubes.server_processes[[length(env$gdalcubes.server_processes) + 1]] <- processx::process$new("gdalcubes_server", c("-b", endpoint, "-p", port[i], "-t", threads), stdin="|")

    Sys.sleep(1.2)
    if (!env$gdalcubes.server_processes[[length(env$gdalcubes.server_processes)]]$is_alive()) {
      stop("gdalcubes_server failed to start")
    }
    message(paste("gdalcubes_server started with pid", env$gdalcubes.server_processes[[length(env$gdalcubes.server_processes)]]$get_pid(), sep=" "))
  }
  invisible()
}

#' Stop all gdalcubes_server processes
#'
#' Stops all gdalcubes_server processes that have been started within the current R session.
#' @export
gcbs_stop_server <- function() {
  stopifnot(requireNamespace("processx", quietly = TRUE))
  x = sum(sapply(.pkgenv$gdalcubes.server_processes, function(x) {
    if (x$is_alive()) {
      x$kill()
      return(1)
    }
    return(0)
    }))
  cat(paste("Stopped", x, "running gdalcubes_server instances\n", sep = " "))
  invisible()
}

#' Status report of gdalcubes_server instances
#'
#' Summarizes the status of all gdalcubes_server processes that have been started in the current R session.
#' @return a data.frame with detailed status information where each row corresponds to one gdalcubes_server process
#' @export
gcbs_server_status <- function() {
  stopifnot(requireNamespace("processx", quietly = TRUE))

  processes <- as.data.frame(t(sapply(.pkgenv$gdalcubes.server_processes, function(p) {
    if (p$is_alive())  {
      return(
        c(name = p$get_name(),
          pid = p$get_pid(),
          cmdline = paste(p$get_cmdline(), collapse =" "),
          alive = p$is_alive(),
          status = p$get_status(),
          start_time = p$get_start_time(),
          cpu_times = as.list(p$get_cpu_times()),
          memory_info_rss = as.numeric(p$get_memory_info()["rss"]),
          memory_info_vmem = as.numeric(p$get_memory_info()["vmem"]),
          username = p$get_username(),
          wd = p$get_wd(),
          exit_status = "-"))
    }
    else {
      return(c(name = "-",
               pid = p$get_pid(),
               cmdline = "-",
               alive = p$is_alive(),
               status = "-",
               start_time = "-",
               cpu_times = rep("-", 4),
               memory_info_rss = "-",
               memory_info_vmem = "-",
               username =  "-",
               wd =  "-",
               exit_status = p$get_exit_status()))
    }
    })))
  return(processes)
}
