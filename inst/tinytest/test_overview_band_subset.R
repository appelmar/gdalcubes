library(gdalcubes)

# A band read from a multi-band image goes through a band-selection VRT whose
# description is empty. When the view is coarse enough to use overviews, the
# source used to be closed before its overview was opened; the overview of ""
# cannot be opened, so the warp went on with a closed dataset and the worker
# crashed. Here: a two-band COG with overviews, mask in band 1 and data in band
# 2, read with a mask at 40 m.

d = tempfile(pattern = "overview_band_subset_")
dir.create(d)
v = cube_view(srs = "EPSG:32609", extent = list(left = 0, right = 6000, bottom = 0, top = 6000,
                                                t0 = "2021-07-05", t1 = "2021-07-05"),
              dx = 10, dy = 10, dt = "P1D")
gdalcubes:::.raster_cube_dummy(v, 2, 1) |>
  apply_pixel(c("band1*4", "band2*500"), names = c("m", "x")) |>
  write_tif(d, prefix = "IMG_", COG = TRUE) |>
  suppressMessages()
fmt = file.path(d, "format.json")
writeLines(c(
  '{"description": "overview test", "tags": ["test"], "pattern": ".*\\\\.tif",',
  ' "images": {"pattern": ".*_([0-9-]{10})\\\\.tif"},',
  ' "datetime": {"pattern": ".*_([0-9-]{10})\\\\.tif", "format": "%Y-%m-%d"},',
  ' "bands": {"SCL": {"pattern": ".*/IMG_.*", "band": 1}, "B04": {"pattern": ".*/IMG_.*", "band": 2}}}'), fmt)
col = create_image_collection(file.path(d, "IMG_2021-07-05.tif"), format = fmt,
                              out_file = file.path(d, "col.db"))

for (p in c(1, 2)) {
  old = gdalcubes_options()$parallel
  gdalcubes_options(parallel = p)
  v40 = cube_view(view = v, dx = 40, dy = 40, resampling = "near")
  out = tempfile(fileext = ".nc")
  r = tryCatch({
    raster_cube(col, v40, mask = image_mask("SCL", values = c(8, 9))) |>
      select_bands("B04") |>
      write_ncdf(out)
    "returned"
  }, error = function(e) conditionMessage(e))
  gdalcubes_options(parallel = old)
  expect_equal(r, "returned", info = paste("masked band-subset read through overviews, parallel", p))
  if (r == "returned") {
    nc = ncdf4::nc_open(out)
    x = ncdf4::ncvar_get(nc, "B04")
    st = ncdf4::ncvar_get(nc, "chunk_status", raw_datavals = TRUE)
    ncdf4::nc_close(nc)
    expect_true(all(st == 0 | st == -2147483647L), info = paste("status OK, parallel", p))
    expect_true(any(!is.na(x)) && all(x[!is.na(x)] == 500), info = paste("data read, parallel", p))
  }
}
