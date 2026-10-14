library(gdalcubes)

# write_ncdf() records a per-chunk status in the output's `chunk_status` variable.
# A read that fails after its image has opened, and a mask image that fails to
# open or to read, must not be recorded as OK: the chunk is otherwise computed
# from the remaining images (or, for a failed mask, from unmasked pixels) with
# nothing to tell the caller.

# Two dates of single-band COGs, a band (B04) and an SCL-style mask. The first
# date is clear (B04 = 500, SCL = 4); the second is cloud (B04 = 900, SCL = 9),
# so with mask values 8 and 9 a value of 900 in the output means that the cloud
# scene was used unmasked.
chunk_status_fixture = function(break_what = "none", srs = "EPSG:32609", size = 400) {
  d = tempfile(pattern = "chunk_status_")
  dir.create(d)
  # COG layout puts the IFD before the tile data, so overwriting a tile corrupts
  # pixels while the file still opens
  tif = function(prefix, date, value) {
    v = cube_view(srs = srs, extent = list(left = 0, right = size, bottom = 0, top = size,
                                           t0 = date, t1 = date),
                  dx = 10, dy = 10, dt = "P1D")
    suppressMessages(write_tif(gdalcubes:::.raster_cube_dummy(v, 1, value), d, prefix = prefix,
                               COG = TRUE, creation_options = list(COMPRESS = "DEFLATE")))
  }
  # write_tif() returns every file in `d` with the prefix, so name them here
  tif("B04_", "2021-07-05", 500); tif("B04_", "2021-07-15", 900)
  tif("SCL_", "2021-07-05", 4); tif("SCL_", "2021-07-15", 9)
  b = file.path(d, c("B04_2021-07-05.tif", "B04_2021-07-15.tif"))
  m = file.path(d, c("SCL_2021-07-05.tif", "SCL_2021-07-15.tif"))
  stopifnot(all(file.exists(c(b, m))))
  fmt = file.path(d, "format.json")
  writeLines(c(
    '{"description": "chunk_status test", "tags": ["test"], "pattern": ".*\\\\.tif",',
    ' "images": {"pattern": ".*_([0-9-]{10})\\\\.tif"},',
    ' "datetime": {"pattern": ".*_([0-9-]{10})\\\\.tif", "format": "%Y-%m-%d"},',
    ' "bands": {"B04": {"pattern": ".*/B04_.*"}, "SCL": {"pattern": ".*/SCL_.*"}}}'), fmt)
  col = create_image_collection(c(b, m), format = fmt, out_file = file.path(d, "col.db"))

  # overwrite the start of the first tile's compressed stream, located from the
  # TileOffsets tag (324) of the first IFD of this little-endian classic TIFF
  corrupt = function(f) {
    con = file(f, "r+b")
    on.exit(close(con))
    stopifnot(identical(readBin(con, "raw", 2), charToRaw("II")))
    seek(con, 4)
    ifd = readBin(con, "integer", size = 4, endian = "little")
    seek(con, ifd)
    n = readBin(con, "integer", size = 2, signed = FALSE, endian = "little")
    tile = NA
    for (i in seq_len(n)) {
      tag = readBin(con, "integer", size = 2, signed = FALSE, endian = "little")
      seek(con, 6, origin = "current")
      val = readBin(con, "integer", size = 4, endian = "little")
      if (tag == 324) tile = val
    }
    stopifnot(!is.na(tile))
    seek(con, tile, rw = "write")
    writeBin(as.raw(rep(0xAB, 64)), con)
  }
  switch(break_what,
         band_read = corrupt(b[2]),
         mask_open = unlink(m[2]),
         mask_read = corrupt(m[2]),
         none = NULL)
  list(col = col, files = c(b, m))
}

# Write a masked B04 cube over the fixture and return its chunk_status (netCDF
# fill, i.e. no status written, read as OK) and its values.
chunk_status_of = function(fx, parallel, srs = "EPSG:32609", extent = c(0, 400, 0, 400),
                           dx = 10, strict = FALSE) {
  old = gdalcubes_options()$parallel
  gdalcubes_options(parallel = parallel)
  on.exit(gdalcubes_options(parallel = old))
  v = cube_view(srs = srs,
                extent = list(left = extent[1], right = extent[2], bottom = extent[3],
                              top = extent[4], t0 = "2021-07-05", t1 = "2021-07-24"),
                dx = dx, dy = dx, dt = "P10D", aggregation = "median", resampling = "near")
  out = tempfile(fileext = ".nc")
  raster_cube(fx$col, v, mask = image_mask("SCL", values = c(8, 9)),
              chunking = c(1, 16, 16), incomplete_ok = !strict) |>
    select_bands("B04") |>
    write_ncdf(out)
  nc = ncdf4::nc_open(out)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  st = as.vector(ncdf4::ncvar_get(nc, "chunk_status", raw_datavals = TRUE))
  st[st == -2147483647L] = 0L
  list(status = st, values = as.vector(ncdf4::ncvar_get(nc, "B04")))
}

for (p in c(1, 2)) {
  # clean reads stay OK, including reprojection, partial coverage and overviews
  r = chunk_status_of(chunk_status_fixture("none"), p)
  expect_true(all(r$status == 0), info = paste("clean read is OK, parallel", p))
  expect_true(all(r$values[!is.na(r$values)] == 500), info = paste("cloud scene masked, parallel", p))
  expect_true(any(!is.na(r$values)), info = paste("clean read has data, parallel", p))

  # the EPSG:32610 fixture at 0..400 lies near x = 668200..668600 in EPSG:32609
  r = chunk_status_of(chunk_status_fixture("none", srs = "EPSG:32610"), p,
                      srs = "EPSG:32609", extent = c(668000, 669000, -300, 700))
  expect_true(all(r$status == 0), info = paste("reprojected read is OK, parallel", p))
  expect_true(any(!is.na(r$values)), info = paste("reprojected read has data, parallel", p))

  r = chunk_status_of(chunk_status_fixture("none"), p, extent = c(-400, 800, -400, 800))
  expect_true(all(r$status == 0), info = paste("partially covered view is OK, parallel", p))

  r = chunk_status_of(chunk_status_fixture("none", size = 6000), p,
                      extent = c(0, 6000, 0, 6000), dx = 40)
  expect_true(all(r$status == 0), info = paste("read through overviews is OK, parallel", p))


  for (b in c("band_read", "mask_open", "mask_read")) {
    fx = chunk_status_fixture(b)
    if (b != "mask_open" && requireNamespace("sf", quietly = TRUE)) {
      # the corrupted file must still open, or this tests an open failure
      info = suppressWarnings(sf::gdal_utils("info", fx$files[if (b == "band_read") 2 else 4],
                                             quiet = TRUE))
      expect_true(length(info) == 1 && !is.na(info) && nzchar(info),
                  info = paste(b, "file still opens"))
    }
    r = suppressWarnings(chunk_status_of(fx, p))
    expect_true(any(r$status != 0), info = paste(b, "is not recorded as OK, parallel", p))
    expect_false(any(r$values == 900, na.rm = TRUE),
                 info = paste(b, "never uses the cloud scene unmasked, parallel", p))
  }

  r = suppressWarnings(chunk_status_of(chunk_status_fixture("mask_open"), p, strict = TRUE))
  expect_true(any(r$status == 1), info = paste("strict mask open failure is ERROR, parallel", p))
}
