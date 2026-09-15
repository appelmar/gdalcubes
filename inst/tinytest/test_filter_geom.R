library(gdalcubes)

# regression test for a segfault when computing a filter_geom() cube
# (https://github.com/appelmar/gdalcubes/issues/110)

v = cube_view(srs = "EPSG:32610",
              extent = list(left = 500000, right = 502000, bottom = 6000000, top = 6002000,
                            t0 = "2020-01-01", t1 = "2020-01-01"),
              dx = 10, dy = 10, dt = "P1D")

wkt = "POLYGON((500500 6000500, 501500 6000500, 501500 6001500, 500500 6001500, 500500 6000500))"

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  filter_geom(wkt, srs = "EPSG:32610") -> cube
expect_equal(dim(cube), c(1, 200, 200))

x = as_array(cube)  # dimensions: band, t, y, x
expect_equal(dim(x), c(1, 1, 200, 200))
expect_equal(sum(!is.na(x)), 10000)     # 100x100 cells inside the polygon
expect_true(all(x[!is.na(x)] == 1))
expect_equal(x[1, 1, 100, 100], 1)      # centre of the polygon
expect_true(is.na(x[1, 1, 10, 10]))     # corners are outside
expect_true(is.na(x[1, 1, 190, 190]))

# with small chunks, chunks completely inside the polygon take a copy fast path;
# the mask must be identical regardless of chunk layout
gdalcubes:::.raster_cube_dummy(v, 1, 1.0, chunking = c(1, 50, 50)) |>
  filter_geom(wkt, srs = "EPSG:32610") -> cube_chunked
y = as_array(cube_chunked)
expect_equal(is.na(x), is.na(y))
