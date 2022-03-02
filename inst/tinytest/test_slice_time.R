library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.05, dy = 0.05)
v

gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  apply_pixel("it", names = "t") |>
  slice_space(c(6.123,49.26933)) |>
  as_array() -> x

expect_equal(as.vector(x), 0:364)  
