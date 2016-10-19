#!/bin/bash


mkdir tmpconda
cd tmpconda
for pkg in latexcodec tabulate monty pybtex palettable spglib pydispatcher pymatgen
do
    conda skeleton pypi $pkg
    conda build $pkg
    anaconda upload /Users/shyuep/miniconda3/conda-bld/osx-64/$pkg-*py35*.tar.bz2
done
cd ..
rm -r tmpconda


mkdir tmpconda
cd tmpconda
for pkg in latexcodec tabulate monty pybtex palettable spglib pydispatcher pymatgen
do
    conda skeleton pypi --python-version 2.7 $pkg
    conda build --python 2.7 $pkg
    anaconda upload --force /Users/shyuep/miniconda3/conda-bld/osx-64/$pkg-*py27*.tar.bz2
done
cd ..
rm -r tmpconda