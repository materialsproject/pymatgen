from setuptools import setup, find_packages

setup (
  name = 'pymatgen',
  version = '1.0.3',
  packages = find_packages(),
  install_requires = ['numpy', 'scipy', 'matplotlib', 'PyCIFRW', 'PyYAML'],

  data_files=[('data', ['pymatgen/core/periodic_table.yaml']),
              ('config', ['pymatgen/io/VaspParameterSets.cfg'])],
  author = 'Shyue Ping Ong, Anubhav Jain, Michael Kocher, Dan Gunter',
  author_email = 'shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, dkgunter@lbl.gov',
  maintainer = 'Shyue Ping Ong',
  summary = 'The Materials Project Python Library',
  url = 'www.materialsproject.org',
  license = 'MIT',
  long_description = 'pymatgen is a Python library for the Materials Project. It includes core structure definition and utilities, electronic structure objects, and convenient IO from VASP and CIF files.',
  keywords = ["vasp", "materials", "project", "electronic", "structure"],
  classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7.2",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Scientists, Developers",
        "License :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Science :: Materials Science",
  ],
  download_url = "http://github.com/downloads/CederGroupMIT/pymatgen_repo/pymatgen_1.0.3.tar.gz",
)

