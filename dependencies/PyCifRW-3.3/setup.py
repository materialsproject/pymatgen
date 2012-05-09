# Setup file for automatic installation of the PyCIFRW
# distribution
from setuptools import setup, Extension

# The compiled scanner for speed

c_scanner = Extension("StarScan",
                      sources = ["lib/lex.yy.c","lib/py_star_scan.c"])

setup(name="PyCifRW",
      version = "3.3",
      description = "CIF/STAR file support for Python",
      author = "James Hester",
      author_email = "jamesrhester at gmail.com",
      url="http://pycifrw.berlios.de",
      py_modules = ['CifFile','yapps3_compiled_rt','YappsStarParser_1_1','YappsStarParser_1_0',
                    'YappsStarParser_DDLm','StarFile'],
      ext_modules = [c_scanner]
      )
