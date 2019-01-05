if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# For details see: https://github.com/rwinlib/gdal2
# For details see: https://github.com/rwinlib/netcdf
if(getRversion() < "3.3.0") setInternet2()

if(!file.exists("../windows/gdal2-2.2.3/include/gdal/gdal.h")){
  download.file("https://github.com/rwinlib/gdal2/archive/v2.2.3.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/netcdf-4.4.1.1-dap/include/netcdf.h")){
  download.file("https://github.com/rwinlib/netcdf/archive/v4.4.1.1-dap.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/libcurl-7.59.0/include/curl/curl.h")){
  download.file("https://github.com/rwinlib/libcurl/archive/v7.59.0.zip", "lib.zip", quiet = TRUE)
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
