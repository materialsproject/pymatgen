FROM python:3.6.2

RUN apt-get update && apt-get install -y gfortran python-openbabel python-vtk python3-tk
RUN curl -sSL https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
 && apt-get update \
 && apt-get install -y git-lfs \
 && rm -rf /var/lib/apt/lists/*
RUN git lfs install