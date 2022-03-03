library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.01, dy = 0.01)
v

# crop space
gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  crop(extent = list(left= 6, right = 9, bottom = 49, top = 52)) -> cube
expect_equal(dim(cube), c(365,300,300))

# crop time
gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  crop(extent = list(t1 = "2021-01-31")) -> cube
expect_equal(dim(cube), c(31,500,500))

gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  crop(extent = list(t0 = "2021-01-30", t1 = "2021-01-31")) -> cube
expect_equal(dim(cube), c(2,500,500))

gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  apply_pixel("it", names = "t") |>
  crop(extent = list(t0 = "2021-01-10", t1 = "2021-01-20")) |>
  as_array() -> x
expect_true(all(x >= 9 & x <= 19))

