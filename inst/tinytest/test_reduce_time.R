library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.02, dy = 0.02)
v


# built-in reducer
gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  reduce_time(c("sum(band1)", "median(band2)", "mean(band3)", "min(band1)","max(band2)", "var(band3)")) |>
  as_array() -> x

expect_true(all(x[1,,,] == 365))
expect_true(all(x[2,,,] == 1))
expect_true(all(x[3,,,] == 1))
expect_true(all(x[4,,,] == 1))
expect_true(all(x[5,,,] == 1))
expect_true(all(x[6,,,] == 0))

gdalcubes:::.raster_cube_dummy(v, 3, 1.0) |>
  reduce_time(c("sum(band1)", "median(band2)"), names=c("A","B")) -> x
expect_true(all(names(x) == c("A","B")))

gdalcubes:::.raster_cube_empty(v, 3) |>
  reduce_time(c("sum(band1)", "median(band2)", "mean(band3)", "min(band1)","max(band2)", "var(band3)")) |>
  as_array() -> x
expect_true(all(is.na(x)))




  

# udf
gdalcubes:::.raster_cube_dummy(v, 2, 1.0) |>
reduce_time(names=c("A", "B"), FUN  = function(x) {
  a = max(x["band1",] + 1:length(x["band1",]))
  b = max(x["band2",] + runif(length(x["band2",])))
  return(c(a, b))
}) |>
as_array() -> x

expect_true(all(x[1,,,] == 366))
expect_true(all(x[2,,,] >= 1 & x[2,,,] <= 2))
