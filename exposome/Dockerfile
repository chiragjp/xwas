## dbmi/exposome-xwas

FROM nampho2/rstudio-server:latest

MAINTAINER "Nam Pho" nam_pho@hms.harvard.edu

RUN yum upgrade -y

RUN yum install libcurl-devel openssl-devel -y

USER rstudio
COPY config /home/rstudio/

RUN mkdir /home/rstudio/.R_libs
RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("nampho2/xwas")'

USER root

#CMD exec /usr/lib/rstudio-server/bin/rserver --server-daemonize 0
