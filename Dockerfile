FROM ubuntu:focal
MAINTAINER Marius Appel <marius.appel@uni-muenster.de>

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y software-properties-common cmake g++ git supervisor wget
RUN apt-get install  -y libudunits2-dev gdal-bin libmagick++-dev pandoc libnetcdf-dev libcurl4-openssl-dev doxygen graphviz libsqlite3-dev libproj-dev libgdal-dev

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt update && apt-get install -y r-base
RUN Rscript -e 'install.packages(c("sf","stars","magick","rmarkdown","ncdf4","Rcpp","jsonlite","RcppProgress", "rstac", "mapview"))'

RUN apt-get install -y gdebi-core
RUN wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.4.1103-amd64.deb
RUN gdebi -n rstudio-server-1.4.1103-amd64.deb

RUN useradd -m -d /home/rstudio rstudio && echo "rstudio:rstudio" | chpasswd

COPY $PWD /opt/gdalcubes_R
RUN cd /opt/gdalcubes_R && R CMD INSTALL .


RUN echo "[supervisord]\nnodaemon=true\nlogfile=/opt/supervisord.log\n[program:rstudio-server]\ncommand=rstudio-server start" > /opt/supervisord.conf


EXPOSE 8787

CMD ["/usr/bin/supervisord", "-c", "/opt/supervisord.conf"]


# docker build -t appelmar/gdalcubes_R .
# run with docker run -d -p 8787:8787 appelmar/gdalcubes_R