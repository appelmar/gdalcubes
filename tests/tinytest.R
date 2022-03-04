
if (getRversion() >= "4.1.0") { # native pipe is needed in tests
  if ( requireNamespace("tinytest", quietly=TRUE) ){
    tinytest::test_package("gdalcubes")
  }
}

