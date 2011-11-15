from setuptools import setup, find_packages

setup (
  name = 'pymatgen',
  version = '1.0.5',
  packages = find_packages(),
  install_requires = ['numpy', 'scipy', 'matplotlib', 'PyCIFRW', 'PyYAML'],

  data_files=[('data', ['pymatgen/core/periodic_table.json']),
              ('config', ['pymatgen/io/VaspParameterSets.cfg'])],
  author = 'Shyue Ping Ong, Anubhav Jain, Michael Kocher, Dan Gunter',
  author_email = 'shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, dkgunter@lbl.gov',
  maintainer = 'Shyue Ping Ong',
  url = 'https://github.com/CederGroupMIT/pymatgen_repo/',
  license = 'MIT',
  long_description = 'pymatgen is a Python library for the Materials Project (www.materialsproject.org). It includes core structure definition and utilities, electronic structure objects, and convenient IO from VASP and CIF files.',
  keywords = ["vasp", "materials", "project", "electronic", "structure"],
  classifiers = [
        "Programming Language :: Python :: 2.7",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
  ],
  download_url = "https://github.com/CederGroupMIT/pymatgen_repo/tarball/master",
)

