# Getting started


## Installation

### Linux source builds
gdalcubes can be compiled from sources via [CMake](https://cmake.org/). CMake automatically checks whether all dependencies (Boost, GDAL, NetCDF, SQLite, and cpprestsdk libraries) are available and 
reports if not. The following commands install gdalcubes from sources. 

```
git clone https://github.com/appelmar/gdalcubes && cd gdalcubes
mkdir -p build 
cd build 
cmake -DCMAKE_BUILD_TYPE=Release ../ 
make 
make install
```

### Windows
All used libraries work under Windows. However, we have not yet tested the compilation on Windows and therefore cannot provide 
detailed instructions or binaries at the moment. You can still use the provided Docker image to run gdalcubes.



### Docker
The `Dockerfile` at the root of the project is built on a minimal Ubuntu installation but installs all dependencies and compiles 
gdalcubes from sources automatically. 


```
git clone https://github.com/appelmar/gdalcubes && cd gdalcubes 
docker build -t appelmar/gdalcubes .
docker run -d -p 11111:1111 appelmar/gdalcubes # runs gdalcubes_server as a deamon 
docker run appelmar/gdalcubes /bin/bash # get a command line where you can run gdalcubes 
``` 



## Sample data

A small (approx. 24 GB ) sample dataset with 20 Landsat images can be downloaded [here](https://uni-muenster.sciebo.de/s/6OjnEyxzt4rk6px/download).






## Next steps

To try out gdalcubes, we recommend reading the [tutorial](tutorial.md).












