library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.01, dy = 0.01)
v

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  apply_pixel("it", names = "t") |>
  slice_time("2021-02-01") |>
  as_array() -> x

expect_true(all(x == 31))  
