from setuptools import setup,find_packages

setup (
  name = 'pymatgen',
  version = '1.0.2',
  packages = find_packages(),

  # Declare your packages' dependencies here, for eg:

  install_requires = ['numpy','scipy','matplotlib','PyCIFRW'],

  author = 'Shyue Ping Ong, Anubhav Jain, Michael Kocher, Dan Gunter',
  author_email = 'shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, dkgunter@lbl.gov',

  summary = 'The Materials Project Python Library',
  url = 'www.materialsproject.org',
  license = '',
  long_description= 'pymatgen is a Python library for the Materials Project. It includes core structure definition and utilities, electronic structure objects, and convenient IO from VASP and CIF files.',

  # could also include long_description, download_url, classifiers, etc.
)
