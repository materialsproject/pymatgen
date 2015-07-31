FROM ubuntu:latest
RUN apt-get -y update && apt-get install -y python-scipy gfortran python-pip python-matplotlib python-openbabel
RUN pip install pymatgen
RUN pip install ipython pyzmq jinja2
RUN pip install tornado jsonschema
CMD python -c 'import pymatgen; print("pymatgen %s installed! Run this image in interactive mode with -t -i to run python or ipython." % pymatgen.__version__)'
