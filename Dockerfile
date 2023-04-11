FROM rocker/rstudio
LABEL maintainer="marius.appel@uni-muenster.de" 

RUN apt update
RUN apt install -y libgdal-dev libcurl4-openssl-dev libnetcdf-dev libudunits2-dev gdal-bin libmagick++-dev pandoc

USER rstudio 
RUN Rscript -e 'install.packages(c("stars","sf","magick","knitr", "rmarkdown", "tinytest", "av","lubridate","gifski", "gdalcubes"))'

USER root
# run e.g. with 
# sudo docker build -t "appelmar/gdalcubes_demo" .
# sudo docker run -d --restart unless-stopped -p 8787:8787 -v $(pwd):/home/rstudio -e PASSWORD=PLEASECHANGE appelmar/gdalcubes_demo
