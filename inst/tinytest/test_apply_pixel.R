library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P6M", 
              dx = 0.01, dy = 0.01)
v

# simple math
gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  apply_pixel(c("band1 - band2", "band1 * 2", "sqrt(band3)", "1", "cos(pi)")) |>
  as_array() -> x

expect_true(all(x[1,,,] == 0))
expect_true(all(x[2,,,] == 2))
expect_true(all(x[3,,,] == 1))
expect_true(all(x[4,,,] == 1))
expect_true(all(x[5,,,] == -1))





# dimension variables
gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  apply_pixel(c("it", "iy", "ix", "left", "right", "top", "bottom", "t0", "t1")) |>
  as_array() -> x

expect_equal(range(x[1,,,]),c(0,1))
expect_equal(range(x[2,,,]),c(0,499))
expect_equal(range(x[3,,,]),c(0,499))
expect_equal(range(x[4,,,]),c(v$space$left,v$space$right - v$space$dx))
expect_equal(range(x[5,,,]),c(v$space$left + v$space$dx,v$space$right))
expect_equal(range(x[6,,,]),c(v$space$bottom + v$space$dy,v$space$top))
expect_equal(range(x[7,,,]),c(v$space$bottom,v$space$top - v$space$dy))
expect_equal(range(as.Date(as.POSIXct(x[8,,,],origin="1970-01-01"))), as.Date(c("2021-01-01", "2021-07-01")))
expect_equal(range(as.Date(as.POSIXct(x[9,,,],origin="1970-01-01"))), as.Date(c("2021-07-01", "2022-01-01")))




