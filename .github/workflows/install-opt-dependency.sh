#!/bin/bash
# Install optional Ubuntu dependencies for pymatgen test workflow

set -e

install_boltztrap() {
    # Install BoltzTraP
    wget -O BoltzTraP.tar.bz2 https://owncloud.tuwien.ac.at/index.php/s/s2d55LYlZnioa3s/download
    tar -jxvf BoltzTraP.tar.bz2

    echo "$(realpath boltztrap-1.2.5/src/)" >> $GITHUB_PATH

    echo "BoltzTraP installation completed."
}

install_vampire() {
    # Install Vampire 5.0
    # TODO: not working, got `FileNotFoundError: [Errno 2] No such file or directory: 'vampire-serial`
    wget https://vampire.york.ac.uk/resources/release-5/vampire-5.0-linux.tar.gz
    tar -zxvf vampire-5.0-linux.tar.gz
    mv linux vampire-5.0

    echo "$(realpath vampire-5.0/)" >> $GITHUB_PATH

    echo "Vampire installation completed."
}

install_babel() {
    # Install openbabel
    # TODO: not working
    sudo apt-get install openbabel
}


# Main script
install_boltztrap
install_vampire
install_babel
