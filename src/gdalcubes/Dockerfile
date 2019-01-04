FROM ubuntu:bionic
MAINTAINER Marius Appel <marius.appel@uni-muenster.de>

RUN apt-get update && apt-get install -y software-properties-common libboost-all-dev cmake g++ libsqlite3-dev git supervisor wget

RUN apt-get install  -y libnetcdf-dev libcurl4-openssl-dev libcpprest-dev doxygen graphviz

# install GDAL from sources, necessary libraries from apt
RUN add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160 && apt-get update && apt-get install -y libproj-dev
RUN apt-get install  -y libxml2-dev libopenjp2-7-dev libhdf4-dev  # install libraries needed for some important drivers

RUN wget https://download.osgeo.org/gdal/2.3.2/gdal-2.3.2.tar.gz && tar -xzf gdal-2.3.2.tar.gz
RUN cd gdal-2.3.2 && ./configure && make -j 2 && make install && ldconfig


# replace with git clone
COPY . /opt/gdalcubes
WORKDIR /opt/gdalcubes
RUN mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release ../ && make -j 2 && make install

COPY supervisord.conf /opt/supervisord.conf

EXPOSE 1111
CMD ["/usr/bin/supervisord", "-c", "/opt/supervisord.conf"]


# docker build -t appelmar/gdalcubes .
# docker run -d -p 11111:1111 appelmar/gdalcubes