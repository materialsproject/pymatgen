FROM ubuntu:latest
ADD "Libraries" "Libraries"
ADD bin bin
RUN apt-get -y update && apt-get install -y python-scipy gfortran python-pip python-matplotlib python-openbabel python-setuptools
RUN export PATH=$PATH:`pwd`/bin
RUN export GULP_LIB=`pwd`/Libraries
RUN pip install pymatgen
RUN pip install ipython pyzmq jinja2
RUN pip install tornado jsonschema
RUN python -c 'import pymatgen; print("pymatgen %s installed! Run this image in interactive mode with -t -i to run python or ipython." % pymatgen.__version__)'
CMD ["ipython"]
