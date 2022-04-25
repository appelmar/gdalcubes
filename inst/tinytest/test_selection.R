library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1M", 
              dx = 0.02, dy = 0.02)


gdalcubes:::.raster_cube_dummy(v, 1, 0,chunking = c(1,67,67)) |>
  apply_pixel(c("ix","iy"), c("ix", "iy")) -> cube
cube[,,240,240] |> as_array() -> x1

gdalcubes:::.raster_cube_dummy(v, 1, 0,chunking = c(1,160,160)) |>
  apply_pixel(c("ix","iy"), c("ix", "iy")) -> cube
cube[,,240,240] |> as_array() -> x2

gdalcubes:::.raster_cube_dummy(v, 1, 0,chunking = c(1,256,256)) |>
  apply_pixel(c("ix","iy"), c("ix", "iy")) -> cube
cube[,,240,240] |> as_array() -> x3

gdalcubes:::.raster_cube_dummy(v, 1, 0,chunking = c(1,500,500)) |>
  apply_pixel(c("ix","iy"), c("ix", "iy")) -> cube
cube[,,240,240] |> as_array() -> x4




expect_true(all(x1[1,,,] == 240))
expect_true(all(x1[2,,,] == 240))

expect_true(all(x2[1,,,] == 240))
expect_true(all(x2[2,,,] == 240))

expect_true(all(x3[1,,,] == 240))
expect_true(all(x3[2,,,] == 240))

expect_true(all(x4[1,,,] == 240))
expect_true(all(x4[2,,,] == 240))


  
gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  crop(extent = list(t0 = "2021-02", t1 = "2021-03")) -> cube
expect_equal(dim(cube), c(2,250,250))

gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  apply_pixel("it", names = "t") |>
  crop(extent = list(t0 = "2021-03", t1 = "2021-10")) |>
  as_array() -> x
expect_true(all(x >= 2 & x <= 9))

