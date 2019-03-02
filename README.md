# gdalcubes_R

[![Build Status](https://travis-ci.org/appelmar/gdalcubes_R.svg?branch=master)](https://travis-ci.org/appelmar/gdalcubes_R)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/appelmar/gdalcubes_R?branch=master&svg=true)](https://ci.appveyor.com/project/appelmar/gdalcubes-r)
[![CRAN](https://www.r-pkg.org/badges/version/gdalcubes)](https://cran.r-project.org/package=gdalcubes)

R package for [gdalcubes](https://github.com/appelmar/gdalcubes): Read and process collections of Earth observation image collection as on demand data cubes.

# Features

- Link to and index existing Earth observation imagery from local files or cloud storage. 
- Read multitemporal, multispectral Earth observation image collections as data cubes by applying on-the-fly reprojection, cropping, and resampling to the desired target cube properties.
- Abstract from complexities in the data like different map projections for adjacent images and different resolutions for different bands.
- Apply R functions on data cube chunks.
- Execute data cube operation chains using parallel processing and lazy evaluation.


# Getting started

- The package includes a vignette that illustrates the basic concepts and functionality on a 
small (< 1 GB) MODIS dataset (see `vignettes/getting_started.Rmd`).

- A tutorial how to use the R package to process Sentinel 2 time series can be found at https://appelmar.github.io/gdalcubes/S2R.html.


# Installation
The package is not (yet) available from CRAN, you have to install it from sources with 

```
library(devtools)
install_git("https://github.com/appelmar/gdalcubes_R", args="--recursive")
```

Please make sure that you are using an up-to-date version of the devtools package and that the [git command line client](https://git-scm.com/downloads) is available on your system. Otherwise, the above command might not clone the gdalcubes C++ library as a submodule under src/gdalcubes.

The package needs a few system libraries ([GDAL](https://www.gdal.org), [NetCDF](https://www.unidata.ucar.edu/software/netcdf), [SQLite](https://www.sqlite.org), [curl](https://curl.haxx.se/libcurl)). 


## Windows

On Windows, you will need [Rtools](https://cran.r-project.org/bin/windows/Rtools). System libraries are automatically downloaded from [rwinlib](https://github.com/rwinlib).


## Linux
Please make sure that you install the system libraries e.g. with the package manager of your Linux distribution. Also make sure that you are using a recent version of GDAL.
On Ubuntu, the following commands wlil install all libraries.

```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libnetcdf-dev libcurl4-openssl-dev libsqlite3-dev libudunits2-dev
```

If the configuration script cannot find the libraries, please install `pkg-config`.



## MacOS
Use [Homebrew](https://brew.sh) to install required libraries with

```
brew install pkg-config
brew install gdal
brew install netcdf
brew install libgit2
brew install udunits
brew install curl
brew install sqlite
```

## Docker
The package includes a Dockerfile that runs [RStudio Server](https://www.rstudio.com/products/rstudio-server/) with the gdalcubes R package. Use the following commands to run a container that becomes accessible from a browser at http://localhost:8787 (user `rstudio`, password `rstudio`).

```
docker build -t appelmar/gdalcubes_R .
docker run -d -p 8787:8787 appelmar/gdalcubes_R
```


# Warning 
The package is still in an early development version. Major changes on the functionality and
names are possible. The documentation is also preliminary and not yet complete.
The package has not been tested on Mac OS.
