if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# For details see: https://github.com/rwinlib/gdal2
# For details see: https://github.com/rwinlib/netcdf

if(!file.exists("../windows/gdal3-3.2.1/include/gdal-3.2.1/gdal.h")){
  download.file("https://github.com/rwinlib/gdal3/archive/v3.2.1.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

# previous winlibs come with sqlite library but not with header file
if(!file.exists("../windows/sqlite-amalgamation-3260000/sqlite3.h")) {
  download.file("https://www.sqlite.org/2018/sqlite-amalgamation-3260000.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
