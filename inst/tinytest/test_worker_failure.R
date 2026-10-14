library(gdalcubes)

# A worker process that dies before writing its log (here: killed by a signal at
# startup) must make write_ncdf() raise, not spin forever reading the missing log
# file. The computation runs in a child R process under a timeout, so a
# regression fails this test instead of hanging the test run.

if (.Platform$OS.type == "unix") {
  d = tempfile(pattern = "worker_failure_")
  dir.create(d)
  worker = file.path(d, "worker_kill.R")
  writeLines(c(
    'j = jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1])',
    'if (j$worker_id == 1) tools::pskill(Sys.getpid(), tools::SIGKILL)',
    'source(system.file("scripts/worker.R", package = "gdalcubes"))'
  ), worker)
  main = file.path(d, "main.R")
  writeLines(c(
    'library(gdalcubes)',
    'ge = asNamespace("gdalcubes")$.pkgenv',
    sprintf('cmd = paste0(file.path(R.home("bin"), "Rscript"), " --no-save --no-restore \\"%s\\"")', worker),
    'gdalcubes:::gc_set_process_execution(2, cmd, FALSE, ge$worker.compression_level,',
    '                                     ge$worker.use_overview_images, ge$worker.gdal_options)',
    'v = cube_view(srs = "EPSG:32609", extent = list(left = 0, right = 4000, bottom = 0, top = 4000,',
    '              t0 = "2021-01-01", t1 = "2021-01-31"), dx = 20, dy = 20, dt = "P1D")',
    'r = tryCatch({ write_ncdf(gdalcubes:::.raster_cube_dummy(v, 1, 1.0), tempfile(fileext = ".nc")); "returned" },',
    '             error = function(e) "raised")',
    'cat("RESULT", r, "\\n")'
  ), main)
  out = suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(main),
                                 stdout = TRUE, stderr = TRUE, timeout = 120))
  expect_true(is.null(attr(out, "status")) || attr(out, "status") != 124,
              info = "write_ncdf() did not hang after a worker was killed")
  expect_true(any(grepl("RESULT raised", out, fixed = TRUE)),
              info = "write_ncdf() raised after a worker was killed")
}
