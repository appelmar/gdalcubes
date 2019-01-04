.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname,pkgname) {

  # call gdalcubes_init()
  if(!Sys.getenv("GDALCUBES_STREAMING") == "1") {
    libgdalcubes_init()
  }
  invisible()
}

.onUnload <- function(libpath) {
  if(!Sys.getenv("GDALCUBES_STREAMING") == "1") {
    libgdalcubes_cleanup()
  }
}


.is_streaming <- function() {
  return(Sys.getenv("GDALCUBES_STREAMING") == "1")
}

.onAttach <- function(libname,pkgname)
{
  # if the package is used inside streaming, redirect stdout to stderr
  # in order to not disturb the (binary) communication with gdalcubes
  if(Sys.getenv("GDALCUBES_STREAMING") == "1") {
    sink(stderr())
  }
  else {
    x = libgdalcubes_version()
    packageStartupMessage(paste("Using gdalcubes library version ", x$VERSION_MAJOR, ".", x$VERSION_MINOR, ".", x$VERSION_PATCH, " (", x$GIT_COMMIT, ")", sep=""))
  }
}
