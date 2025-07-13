#!/usr/bin/env bash
set -euo pipefail
set -x

# Use a custom bin directory for easier GitHub Actions caching
BIN_DIR="${1:-$PWD/opt/bin}" && mkdir -p "$BIN_DIR"

# Install BoltzTraP
wget --no-verbose -O BoltzTraP.tar.bz2 https://owncloud.tuwien.ac.at/index.php/s/s2d55LYlZnioa3s/download
tar -jxf BoltzTraP.tar.bz2
cp boltztrap-1.2.5/src/{BoltzTraP,x_trans} "$BIN_DIR"

# Install Vampire 5.0
# TODO: Accepts self-signed cert (https://github.com/richard-evans/vampire/issues/122)
wget --no-verbose --no-check-certificate https://vampire.york.ac.uk/resources/release-5/vampire-5.0-linux.tar.gz
tar -zxf vampire-5.0-linux.tar.gz && mv linux vampire-5.0
cp vampire-5.0/vampire-serial "$BIN_DIR"

# Install Voro++ and ZEO++
wget --no-verbose http://www.zeoplusplus.org/zeo++-0.3.tar.gz
tar -xzf zeo++-0.3.tar.gz
make -C zeo++-0.3/voro++/src -s CFLAGS="-w"
make -C zeo++-0.3 -s CFLAGS="-w"
cp zeo++-0.3/voro++/src/voro++ "$BIN_DIR"
cp zeo++-0.3/network "$BIN_DIR"

# Install mcsqs (from ATAT)
# TODO: cannot build on MacOS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  wget --no-verbose https://axelvandewalle.github.io/www-avdw/atat/atat3_50.tar.gz
  tar -zxf atat3_50.tar.gz

  # Replace `BINDIR` in makefile
  sed -i "s|^BINDIR=.*|BINDIR=$BIN_DIR|" atat/makefile

  sudo apt install csh -y -q
  make -C atat && make -C atat install
fi

# Add to path
echo "$BIN_DIR" >> "$GITHUB_PATH"
