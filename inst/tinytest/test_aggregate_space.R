library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1M", 
              dx = 1, dy = 1)
v


gdalcubes:::.raster_cube_dummy(v, 1, 1.0) |>
  aggregate_space(dx = 2, dy = 2, method = "sum") |>
  as_array() -> x

expect_true(all(x[1,1,,] == matrix(c(1,2,2,2,4,4,2,4,4), nrow = 3, byrow = TRUE)))
expect_true(all(x[1,12,,] == matrix(c(1,2,2,2,4,4,2,4,4), nrow = 3, byrow = TRUE)))



gdalcubes:::.raster_cube_dummy(v, 1, 1.0, chunking = c(1,4,4)) |>
  aggregate_space(dx = 2, dy = 2, method = "count") |>
  as_array() -> y
expect_equal(x,y)


gdalcubes:::.raster_cube_dummy(v, 1, 1.0, chunking = c(1,2,10)) |>
  aggregate_space(dx = 2, dy = 2, method = "count") |>
  as_array() -> z
expect_equal(x,z)


gdalcubes:::.raster_cube_dummy(v, 3, 2.0, chunking = c(1,2,10)) |>
  aggregate_space(dx = 3, dy = 3, method = "mean") |>
  as_array() -> z
expect_true(all(z == 2))


gdalcubes:::.raster_cube_dummy(v, 3, 2.0, chunking = c(3,4,4)) |>
  aggregate_space(dx = 3, dy = 3, method = "mean") |>
  as_array() -> z
expect_true(all(z == 2))


gdalcubes:::.raster_cube_dummy(v, 3, 2.0, chunking = c(3,4,4)) |>
  aggregate_space(dx = 3, dy = 3, method = "var") |>
  as_array() -> z
expect_true(all(z == 0))


gdalcubes:::.raster_cube_dummy(v, 3, 2.0, chunking = c(1,3,3)) |>
  aggregate_space(dx = 3, dy = 3, method = "sd") |>
  as_array() -> z
expect_true(all(z == 0))



### 
# library(terra)
# x <- rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0,ymax=10)
# values(x) <- m <- matrix(1:100, nrow=10, ncol=10, byrow = TRUE)
# writeRaster(x, "m.tif",overwrite=TRUE)
# plot(x)


#xx |> aggregate_space(dx=2.7, dy=2.7, method = "median") |> 
a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,10,10)) |> 
  aggregate_space(fact=5, method = "mean") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,28,73,78), nrow=2, byrow=TRUE))
  

a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,3,3)) |> 
  aggregate_space(fact=5, method = "mean") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,28,73,78), nrow=2, byrow=TRUE))

a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,1,1)) |> 
  aggregate_space(fact=5, method = "mean") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,28,73,78), nrow=2, byrow=TRUE))


a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,3,3)) |> 
  aggregate_space(dx=4.2, dy = 4.2, method = "max") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,27,30, 63, 67, 70, 93, 97, 100), nrow=3, byrow=TRUE))


a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,9,2)) |> 
  aggregate_space(dx=4.2, dy = 4.2, method = "max") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,27,30, 63, 67, 70, 93, 97, 100), nrow=3, byrow=TRUE))

a = stack_cube( "m.tif",bands = "B1",datetime_values = "2000-01-01",chunking = c(1,1,1)) |> 
  aggregate_space(dx=4.2, dy = 4.2, method = "max") |> 
  as_array() |> drop()
expect_equivalent(a, matrix(c(23,27,30, 63, 67, 70, 93, 97, 100), nrow=3, byrow=TRUE))


