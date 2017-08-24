#ubuntu:16.04
FROM continuumio/anaconda
#RUN apt-get update && apt-get install -y wget make git python default-jre default-jdk python-dev python-matplotlib gcc swig python-numpy libgsl2 \
#    gsl-bin pkg-config zlib1g-dev libcairo2-dev libnetpbm10-dev netpbm libpng12-dev libjpeg-dev python-pyfits zlib1g-dev \
#    libbz2-dev libcfitsio3-dev python-photutils python-pip bedtools bbtools
RUN conda install -c anaconda graphviz && conda install -c bioconda bedtools bbmap
#RUN wget https://anaconda.org/jlevy44/subgenomeClustering/2017.08.24.0831/download/subgenomeClustering.yml && conda env create -f subgenomeClustering.yml
#CMD ["source", "activate", "subgenomeClustering"]
#RUN conda install -c anaconda graphviz
