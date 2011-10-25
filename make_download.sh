#!/bin/bash

version=`awk '/version/ {print $3}' setup.py`
version=${version:1}
version=${version%\',}

dir=pymatgen_$version
mkdir $dir

cp -r pymatgen $dir
cp README.md $dir
cp run_me_first.sh $dir
cp setup.py $dir
find $dir -name "tests" -exec rm -r '{}' \;
find $dir -name "pymatgen.cfg" -exec rm -r '{}' \;
find $dir -name "*.pyc" -exec rm -r '{}' \;
find $dir -name ".DS_Store" -exec rm -r '{}' \;

tar -zcf $dir.tar.gz $dir
rm -rf $dir
