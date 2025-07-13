#!/usr/bin/env bash
set -euo pipefail
set -x

# Install BoltzTraP
wget -O BoltzTraP.tar.bz2 https://owncloud.tuwien.ac.at/index.php/s/s2d55LYlZnioa3s/download
tar -jxf BoltzTraP.tar.bz2
echo "$(realpath boltztrap-1.2.5/src)" >> "$GITHUB_PATH"

# Install Vampire 5.0
wget https://vampire.york.ac.uk/resources/release-5/vampire-5.0-linux.tar.gz
tar -zxf vampire-5.0-linux.tar.gz
mv linux vampire-5.0
echo "$(realpath vampire-5.0)" >> "$GITHUB_PATH"

# Install Voro++ and ZEO++
wget http://www.zeoplusplus.org/zeo++-0.3.tar.gz
tar -xzf zeo++-0.3.tar.gz

make -C zeo++-0.3/voro++/src -s CFLAGS="-w"
echo "$(realpath zeo++-0.3/voro++/src)" >> "$GITHUB_PATH"

make -C zeo++-0.3 -s CFLAGS="-w"
echo "$(realpath zeo++-0.3)" >> "$GITHUB_PATH"

# TODO: Install mcsqs (from ATAT)
# wget https://axelvandewalle.github.io/www-avdw/atat/atat3_50.tar.gz
# tar -zxf atat3_50.tar.gz && cd atat
