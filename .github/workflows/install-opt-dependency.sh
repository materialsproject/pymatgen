#!/bin/bash
# Install optional Ubuntu dependencies for pymatgen test workflow

set -e

install_boltztrap() {
    # Install BoltzTraP
    wget -O BoltzTraP.tar.bz2 https://owncloud.tuwien.ac.at/index.php/s/s2d55LYlZnioa3s/download
    tar -jxvf BoltzTraP.tar.bz2

    ln -s boltztrap-1.2.5/src/x_trans /usr/local/bin/
    ln -s boltztrap-1.2.5/src/BoltzTraP /usr/local/bin/

    echo "BoltzTraP installation completed."
}

# Main script
install_boltztrap
