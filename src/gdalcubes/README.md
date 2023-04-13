# gdalcubes - Earth observation data cubes from GDAL image collections 

[![Build Status](https://travis-ci.org/appelmar/gdalcubes.svg?branch=master)](https://travis-ci.org/appelmar/gdalcubes)

<a href="https://gdalcubes.github.io/docs"><img src="gdalcubes_logo_small.png" align="right" height="180"/></a>
**gdalcubes** is a library to represent collections of Earth Observation (EO) images as _on demand_ data cubes (or _multidimensional arrays_). Users define data cubes by spatiotemporal extent, resolution, and spatial reference system and let gdalcubes read only relevant parts of the data and simultaneously apply reprojection, resampling, and cropping (using [gdalwarp](https://www.gdal.org/gdalwarp.html)). Data cubes may be simply exported as NetCDF files or directly streamed chunk-wise into external software such as R or Python. The library furthermore implements simple operations to reduce data cubes over time, to apply pixel-wise arithmetic expressions, and to filter by space, time, and bands. 

gdalcubes is not a database, i.e., it does not need to store additional copies of the imagery but instead
simply links to and indexes existing files / GDAL datasets. Using [GDAL virtual file systems](https://www.gdal.org/gdal_virtual_file_systems.html), it can directly access
data in cloud storage and run computations in distributed environments with `gdalcubes_server` and Docker. 

The library is written in C++ and includes a basic command line interface and an [R package]((https://github.com/appelmar/gdalcubes_R)). A python package is
planned for the future. gdalcubes is licensed under the [MIT license](https://opensource.org/licenses/MIT).

# Core features:

- Create image collections that link to and index existing imagery from local files or cloud storage 
- Read multitemporal, multispectral image collections as on demand data cubes with desired spatiotemporal resolution, extent, and map projection
- Abstract from complexities in the data like different map projections for adjacent images and different resolutions for different bands
- Stream chunks of data cubes to external programs (e.g. R, python)
- Scale computations on data cubes in distributed environments with `gdalcubes_server` and Docker (yet experimental)



# Installation


## Installation from sources

gdalcubes uses CMake and can be compiled with a typical CMake workflow as listed below.

```
git clone https://github.com/appelmar/gdalcubes && cd gdalcubes
mkdir -p build 
cd build 
cmake -DCMAKE_BUILD_TYPE=Release ../ -DCMAKE_INSTALL_PREFIX=/usr
make 
sudo make install
```

You might need to install a few libraries before compiling gdalcubes successfully. On Ubuntu `apt install libgdal-dev libnetcdf-dev libcurl4-openssl-dev libsqlite3-dev` will install all libraries needed to compile 
the core gdalcubes library. If you want to compile the command line interface, you will furthermore need `apt install libboost-program-options-dev libboost-system-dev`
and running gdalcubes as a server additionally requires `apt install libcpprest-dev`.





## Docker image
This repository includes a Docker image which you can use either to run the gdalcubes command line interface interactively
or to run `gdalcubes_server` as a service for distributed processing. The commands below demonstrate how to build the image and run a container.
Notice that the image builds GDAL from sources, which might take up to 30 minutes.
 

```
git clone https://github.com/appelmar/gdalcubes && cd gdalcubes 
docker build -t appelmar/gdalcubes .
docker run -d -p 11111:1111 appelmar/gdalcubes # runs gdalcubes_server as a deamon 
docker run appelmar/gdalcubes /bin/bash # get a command line where you can run gdalcubes 
``` 


# R package
The [gdalcubes R package](https://github.com/appelmar/gdalcubes_R) is hosted on https://github.com/appelmar/gdalcubes_R.
It includes a Dockerfile that runs [RStudio Server](https://www.rstudio.com/products/rstudio-server/) with the gdalcubes R package.

# Documentation
More detailed information can be found at the documentation page under https://gdalcubes.github.io/. 

# Warning
The library is still in an early development version. Major changes are possible to make gdalcubes more user-friendly, more stable, faster, and more robust.
The documentation is also preliminary and not yet complete.

# Credits
gdalcubes uses the following open source libraries. Detailed licensing and copyright information can be found at https://gdalcubes.github.io/source/introduction/credits.html and in LICENSE_THIRDPARTY.


**[GDAL](https://www.gdal.org/):  A translator library for raster and vector geospatial data formats**

**[json11](https://github.com/dropbox/json11)**

**[SQLite](https://www.sqlite.org/): A self-contained, high-reliability, embedded, full-featured, public-domain, SQL database engine**

**[CURL](https://curl.haxx.se/): Command line tool and library for transferring data with URLs**

**[TinyExpr](https://github.com/codeplea/tinyexpr): A very small recursive descent parser and evaluation engine for math expressions**

**[netCDF](https://www.unidata.ucar.edu/software/netcdf): The Unidata network Common Data Form C library**
   
**[tiny-process-library](https://gitlab.com/eidheim/tiny-process-library): A small platform independent library making it simple to create and stop new processes in C++**

**[Catch2](https://github.com/catchorg/Catch2): A modern, C++-native, header-only, test framework for unit-tests, TDD and BDD**
       
**[Date](https://github.com/HowardHinnant/date): A date and time library based on the C++11/14/17 <chrono> header**   

**[cpprestsdk](https://github.com/Microsoft/cpprestsdk)**

**[Boost.Filesystem](https://www.boost.org/doc/libs/1_68_0/libs/filesystem/doc/index.htm)**

**[Boost.Program_options](https://www.boost.org/doc/libs/1_68_0/doc/html/program_options.html)**
