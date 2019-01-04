# gdalcubes - Earth observation data cubes from GDAL image collections

**gdalcubes** is a library to represent collections of Earth Observation (EO) images
as _on demand_ data cubes (or _multidimensional arrays_). Users define data cubes by spatiotemporal extent, resolution, and 
spatial reference system and let gdalcubes read only relevant parts of the data and simultaneously apply reprojection, resampling, and cropping (using [gdalwarp](https://www.gdal.org/gdalwarp.html)).
Data cubes may be simply exported as NetCDF files or directly streamed chunk-wise into external software such as R or Python. The library furthermore
implements simple operations to reduce data cubes over time, to apply pixel-wise arithmetic expressions, and to filter by space, time, and bands.

gdalcubes is not a database, i.e., it does not need to store additional copies of the imagery but instead
simply links to and indexes existing files / GDAL datasets. Using [GDAL virtual file systems](https://www.gdal.org/gdal_virtual_file_systems.html), it can directly access
data in cloud storage and run computations in distributed environments with `gdalcubes_server`. 

The library is written in C++ and includes a basic command line interface as well as a package for R. A python package is
planned for the future. gdalcubes is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

## Core features:

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
make install
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

The Dockerfile at `demo/Dockerfile` furthermore runs [RStudio Server](https://www.rstudio.com/products/rstudio-server/) with the gdalcubes R package.


```
# Make sure that you have successfully built the docker image above before
docker build -t appelmar/gdalcubes_demo .
docker run -d -p 8787:8787 appelmar/gdalcubes_demo 

# -> open http://localhost:8787 in a browser and login with user/password rstudio/rstudio.
``` 




## R package
You can install the current R package version directly from GitHub using the `devtools` package as shown below.
Notice that this works only if you have successfully installed the gdalcubes library on your system before.

```
library(devtools)
install_github("appelmar/gdalcubes", subdir="Rpkg/gdalcubes")
```

The package includes a vignette that illustrates the basic ideas on a small MODIS dataset.


# Documentation
More detailed information can be found at the documentation page under https://appelmar.github.io/gdalcubes. 



# Credits

gdalcubes uses the following great open source libraries. Detailed licensing and copyright information can be found at https://appelmar.github.io/gdalcubes/credits.html.


**[GDAL](https://www.gdal.org/):  A translator library for raster and vector geospatial data formats**

**[json](https://github.com/nlohmann/json): JSON for Modern C++**

**[SQLite](https://www.sqlite.org/): A self-contained, high-reliability, embedded, full-featured, public-domain, SQL database engine**

**[CURL](https://curl.haxx.se/): Command line tool and library for transferring data with URLs**

**[ExprTk](http://www.partow.net/programming/exprtk/): A C++ Mathematical Expression Parsing and Evaluation Library**
 
**[netCDF](https://www.unidata.ucar.edu/software/netcdf): The Unidata network Common Data Form C library**
   
**[tiny-process-library](https://gitlab.com/eidheim/tiny-process-library): A small platform independent library making it simple to create and stop new processes in C++**

**[Catch](https://www.gdal.org/): A modern, C++-native, header-only, test framework for unit-tests, TDD and BDD**
       
**[Date](https://github.com/HowardHinnant/date): A date and time library based on the C++11/14/17 <chrono> header**   

**[cpprestsdk](https://github.com/Microsoft/cpprestsdk)**

**[Boost.Filesystem](https://www.boost.org/doc/libs/1_68_0/libs/filesystem/doc/index.htm)**

**[Boost.Program_options](https://www.boost.org/doc/libs/1_68_0/doc/html/program_options.html)**
