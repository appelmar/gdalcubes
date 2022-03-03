library(gdalcubes)
v = cube_view(srs = "EPSG:4326", extent = list(left = 5, right = 10, bottom = 48, top = 53, 
                                               t0 = "2021-01-01", t1 = "2021-12-31"), dt = "P1D", 
              dx = 0.01, dy = 0.01)

expect_equal(v$space$nx, 500)
expect_equal(v$space$ny, 500)
expect_equal(v$space$dx, 0.01)
expect_equal(v$space$dy, 0.01)
expect_equal(v$time$nt, 365)

v = cube_view(v, dt = "P1M")
expect_equal(v$space$nx, 500)
expect_equal(v$space$ny, 500)
expect_equal(v$space$dx, 0.01)
expect_equal(v$space$dy, 0.01)
expect_equal(v$time$nt, 12)
expect_equal(v$time$t0, "2021-01")
expect_equal(v$time$t1, "2021-12")

v = cube_view(v, nx = 5000, ny = 5000)
expect_equal(v$space$nx, 5000)
expect_equal(v$space$ny, 5000)
expect_equal(v$space$dx, 0.001)
expect_equal(v$space$dy, 0.001)
