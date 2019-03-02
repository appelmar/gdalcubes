.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname,pkgname) {

  # call gdalcubes_init()
  if(!Sys.getenv("GDALCUBES_STREAMING") == "1") {
    libgdalcubes_init()
    libgdalcubes_add_format_dir(file.path(system.file(package="gdalcubes"),"formats")) # add collection formats directory 
  }
  
  .pkgenv$compression_level = 0
  .pkgenv$cube_cache = new.env()
  .pkgenv$use_cube_cache = TRUE
  
  # for windows, rwinlib includes GDAL data and PROJ data in the package and we must set the environment variables
  # PROJ_LIB and GDAL_DATA to make sure GDAL finds the data at the package location
  if (dir.exists(system.file("proj", package="gdalcubes"))) {
    Sys.setenv("PROJ_LIB" = system.file("proj", package="gdalcubes"))
  }
  if (dir.exists(system.file("gdal", package="gdalcubes"))) {
    Sys.setenv("GDAL_DATA" = system.file("gdal", package="gdalcubes"))
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
    packageStartupMessage(paste("Using gdalcubes library version ", x$VERSION_MAJOR, ".", x$VERSION_MINOR, ".", x$VERSION_PATCH, sep=""))
  }
}
