library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 15, bottom = 48, top = 58, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P365D", 
              dx = 1, dy = 1)
v

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  window_space("max(band1)","min(band1)", window = c(5,5)) |>
  as_array() -> x

expect_true(all(x == 1))  



gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  window_space("count(band1)", window = c(3,3)) |>
  as_array() -> x
Xtrue = matrix(9, 10,10)
Xtrue[c(1,10),] = Xtrue[,c(1,10)] = 6
Xtrue[1,1] = Xtrue[1,10] = Xtrue[10,1] = Xtrue[10,10] = 4
expect_true(all(x[1,1,,] == Xtrue))


K = matrix(1, 3,3)

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  window_space(kernel = K, pad = 0) |>
  as_array() -> x
expect_true(all(x[1,1,,] == Xtrue))

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  window_space(kernel = K, pad = "REFLECT") |>
  as_array() -> x
expect_true(all(x == 9))


gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  window_space(kernel = K, pad = "REPLICATE") |>
  as_array() -> x
expect_true(all(x == 9))




gdalcubes:::.raster_cube_dummy(v, 1, 1.0, chunking = c(1,3,2)) |>
  window_space(kernel = K, pad = "REPLICATE") |>
  as_array() -> x
expect_true(all(x == 9))

