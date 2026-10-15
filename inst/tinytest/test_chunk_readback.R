library(gdalcubes)

# The main process reads each worker's chunk file back with chunk_data::read_ncdf().
# A file it cannot read must make the chunk an ERROR, not an OK chunk of NA (or of
# an unread buffer merged as data). Worker 1 is replaced by a script that writes
# a broken file for each of its chunks and exits 0; worker 0 runs normally, so
# its chunks are the clean control.
#
# Not covered: a failing nc_get_var_double(). Workers write netCDF-4; a truncated
# netCDF-4 file already fails at nc_open(), and one with damaged data bytes reads
# back without an error, so neither reaches that check from here.

if (.Platform$OS.type != "unix") exit_file("stand-in worker script needs Unix")

d = tempfile(pattern = "chunk_readback_")
dir.create(d)
worker = file.path(d, "worker_broken.R")
writeLines(c(
  'j = jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1])',
  '# a worker only writes its log when it logs something, and on a worker failure the main',
  '# process waits forever for a missing log (fixed separately), so create it up front',
  'invisible(file.create(j$gdalcubes_options$log_file))',
  'if (j$worker_id == 0) source(system.file("scripts/worker.R", package = "gdalcubes"))',
  'if (j$worker_id == 1) {',
  '  kind = Sys.getenv("GDALCUBES_TEST_READBACK")',
  '  for (id in seq(1, as.integer(Sys.getenv("GDALCUBES_TEST_NCHUNKS")) - 1, by = 2)) {',
  '    tmp = file.path(j$workdir, paste0(".broken_", id, ".nc"))',
  '    if (kind == "not_netcdf") writeLines("not a netCDF file", tmp)',
  '    else {',
  '      # a real worker file carries chunk_status = 0 and dimensions b, t, y, x',
  '      dn = if (kind == "no_dimensions") c("u", "v", "w", "z") else c("x", "y", "t", "b")',
  '      dims = mapply(function(n, len) ncdf4::ncdim_def(n, "", seq_len(len), create_dimvar = FALSE),',
  '                    dn, c(20, 20, 1, 1), SIMPLIFY = FALSE)',
  '      vn = if (kind == "no_variable") "other" else "value"',
  '      nc = ncdf4::nc_create(tmp, list(ncdf4::ncvar_def(vn, "", dims, NA, prec = "double")))',
  '      ncdf4::ncatt_put(nc, 0, "chunk_status", 0L, prec = "int")',
  '      ncdf4::ncvar_put(nc, vn, rep(900, 400))',
  '      ncdf4::nc_close(nc)',
  '    }',
  '    file.rename(tmp, file.path(j$workdir, paste0(id, ".nc")))',
  '  }',
  '}'
), worker)

v = cube_view(srs = "EPSG:32609", extent = list(left = 0, right = 400, bottom = 0, top = 400,
                                                t0 = "2021-07-05", t1 = "2021-07-05"),
              dx = 10, dy = 10, dt = "P1D")
nchunks = 4  # 40 x 40 cells in 20 x 20 chunks: worker 0 computes 0 and 2, worker 1 "computes" 1 and 3

readback_status = function(kind) {
  Sys.setenv(GDALCUBES_TEST_READBACK = kind, GDALCUBES_TEST_NCHUNKS = nchunks)
  ge = gdalcubes:::.pkgenv
  cmd = paste0(file.path(R.home("bin"), "Rscript"), " --no-save --no-restore \"", worker, "\"")
  gdalcubes:::gc_set_process_execution(2, cmd, FALSE, ge$worker.compression_level,
                                       ge$worker.use_overview_images, ge$worker.gdal_options)
  on.exit({
    gdalcubes_options(parallel = gdalcubes_options()$parallel)  # restore the default worker
    Sys.unsetenv(c("GDALCUBES_TEST_READBACK", "GDALCUBES_TEST_NCHUNKS"))
  })
  out = tempfile(fileext = ".nc")
  write_ncdf(gdalcubes:::.raster_cube_dummy(v, 1, 500, chunking = c(1, 20, 20)), out)
  nc = ncdf4::nc_open(out)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  st = as.vector(ncdf4::ncvar_get(nc, "chunk_status", raw_datavals = TRUE))
  st[st == -2147483647L] = 0L  # netCDF fill: no status written, read as OK
  list(status = st, values = as.vector(ncdf4::ncvar_get(nc, "band1")))
}

for (kind in c("not_netcdf", "no_dimensions", "no_variable")) {
  r = readback_status(kind)
  expect_equal(r$status[c(2, 4)], c(1L, 1L), info = paste(kind, "chunks are recorded as ERROR"))
  expect_equal(r$status[c(1, 3)], c(0L, 0L), info = paste(kind, "clean chunks stay OK"))
  expect_equal(sum(r$values == 500, na.rm = TRUE), 800L, info = paste(kind, "clean chunks hold their data"))
}
