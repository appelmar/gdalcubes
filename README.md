# gdalcubes_R
R package for [gdalcubes](https://github.com/appelmar/gdalcubes) to process collections of Earth observation image collection as on demand data cubes.

# Features

- Link to and index existing Earth observation imagery from local files or cloud storage. 
- Read multitemporal, multispectral Earth observation image collections as data cubes by applying on-the-fly reprojection, cropping, and resampling to the desired target cube properties.
- Abstract from complexities in the data like different map projections for adjacent images and different resolutions for different bands.
- Apply R functions on data cube chunks.
- Execute data cube operation chains using parallel processing and lazy evaluation.


# Installation
The package is not available from CRAN, you have to install it from sources. On Windows, you will need [Rtools](https://cran.r-project.org/bin/windows/Rtools). Package building automatically downloads dependencies from [rwinlib](https://github.com/rwinlib) on Windows. For Linux, you will need the following system libraries:

- GDAL
- NetCDF 
- SQLite  
- CURL 

Since the package includes the gdalcubes C++ library as a git submodule in its `src` folder, you must install the [git command line client](https://git-scm.com/downloads) before you can use the devtools package for installation as shown below.

```
library(devtools)
install_git("https://github.com/appelmar/gdalcubes_R", args="--recursive")
```

# Docker
The package includes a Dockerfile that runs [RStudio Server](https://www.rstudio.com/products/rstudio-server/) with the gdalcubes R package. Use the following commands to run a container that becomes accessible from a browser at http://localhost:8787 (user `rstudio`, password `rstudio`).

```
docker build -t appelmar/gdalcubes_R .
docker run -d -p 8787:8787 appelmar/gdalcubes_R
```

# Getting started
The package includes a vignette that illustrates the basic concepts and functionality on a 
small (< 1 GB) MODIS dataset (see `vignettes/getting_started.Rmd`).


# Warning 
The package is still in an early development version. Major changes on the functionality and
names are possible. The documentation is also preliminary and not yet complete.
The package has not been tested on Mac OS.
