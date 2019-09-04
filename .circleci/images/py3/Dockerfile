FROM python:3.7.3

RUN apt-get update && apt-get install -y gfortran
# python-openbabel python-vtk python3-tk
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
RUN export PATH=$HOME/miniconda3/bin:$PATH