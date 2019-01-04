FROM appelmar/gdalcubes

# install R and a few packages
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt update -q && DEBIAN_FRONTEND=noninteractive apt install -q -y r-base
RUN Rscript -e "install.packages(c('devtools'))"

RUN apt install -y gdebi-core
RUN wget https://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb
RUN gdebi -n rstudio-server-1.1.463-amd64.deb

RUN useradd -m -d /home/rstudio rstudio && echo "rstudio:rstudio" | chpasswd

# TODO: replace by installing from github

RUN Rscript -e "install.packages(c('RcppProgress', 'ncdf4'))"
RUN cd /opt/gdalcubes/Rpkg &&  R CMD INSTALL gdalcubes

RUN echo "[program:rstudio-server]\ncommand=rstudio-server start" >> /opt/supervisord.conf


EXPOSE 8787

#CMD /bin/bash
CMD ["/usr/bin/supervisord", "-c", "/opt/supervisord.conf"]


# docker build -t appelmar/gdalcubes_demo .
# run with docker run -d -p 8787:8787 -p 1111:1111  appelmar/gdalcubes_demo