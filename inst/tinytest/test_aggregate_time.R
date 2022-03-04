library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.02, dy = 0.02)
v


gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  aggregate_time(dt = "P1M", method = "sum") |>
  as_array() -> x
expect_equal(dim(x), c(1, 12, 250, 250))

expect_true(all(x[1,1,,] == 31))
expect_true(all(x[1,2,,] == 28))
expect_true(all(x[1,3,,] == 31))
expect_true(all(x[1,4,,] == 30))
expect_true(all(x[1,5,,] == 31))
expect_true(all(x[1,6,,] == 30))
expect_true(all(x[1,7,,] == 31))
expect_true(all(x[1,8,,] == 31))
expect_true(all(x[1,9,,] == 30))
expect_true(all(x[1,10,,] == 31))
expect_true(all(x[1,11,,] == 30))
expect_true(all(x[1,12,,] == 31))


gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  aggregate_time(dt = "P1M", method = "count") |>
  as_array() -> y
expect_equal(x,y)



gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  aggregate_time(dt = "P5D", method = "mean") |>
  as_array() -> x
expect_equal(dim(x), c(1,73,250,250))
expect_true(all(x == 1))


gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  aggregate_time(fact = 5, method = "median") |>
  as_array() -> y
expect_equal(dim(x),dim(y))
expect_true(all(y == 1))

# TODO: count reducer might lead to NA instead of 0 on
# some platforms, avoid test for now
# gdalcubes:::.raster_cube_empty(v, 1, 1.0) |>
#   aggregate_time(dt = "P1M", method = "count") |>
#   as_array() -> x
# expect_true(all(x == 0))


gdalcubes:::.raster_cube_empty(v, 1, 1.0) |>
  aggregate_time(dt = "P1M", method = "sum") |>
  as_array() -> x
expect_true(all(is.na(x)))



