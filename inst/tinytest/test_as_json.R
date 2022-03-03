library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.01, dy = 0.01)
v



gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  reduce_space(c("sum(band1)", "median(band2)", "mean(band3)", "min(band1)","max(band2)", "var(band3)")) -> cube
cube |>as_json() |> json_cube() -> cube1

expect_equal(names(cube), names(cube1))
expect_equal(dimensions(cube), dimensions(cube1))


# TODO: add for all cube types