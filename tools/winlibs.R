if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# For details see: https://github.com/rwinlib/gdal3

if(!file.exists("../windows/gdal3-3.4.1/include/gdal.h")){
  download.file("https://github.com/rwinlib/gdal3/archive/v3.4.1.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
