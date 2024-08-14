#!/bin/bash
# Install optional Ubuntu dependencies for pymatgen test workflow

set -e

install_boltztrap() {
    # Install BoltzTraP
    wget -q -O BoltzTraP.tar.bz2 https://owncloud.tuwien.ac.at/index.php/s/s2d55LYlZnioa3s/download
    tar -jxf BoltzTraP.tar.bz2

    echo "$(realpath boltztrap-1.2.5/src/)" >> $GITHUB_PATH

    echo "BoltzTraP installation completed."
}

install_vampire() {
    # Install Vampire 5.0
    wget -q https://vampire.york.ac.uk/resources/release-5/vampire-5.0-linux.tar.gz
    tar -zxf vampire-5.0-linux.tar.gz
    mv linux vampire-5.0

    echo "$(realpath vampire-5.0/)" >> $GITHUB_PATH

    echo "Vampire installation completed."
}


# Main script
install_boltztrap
install_vampire
