# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:/global/home/groups/lr_mp/ewcspottesmith/bin
PATH=$PATH:/clusterfs/mp/software/cp2k/plumed2.0_install/bin

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/clusterfs/mp/software/cp2k/plumed2.0_install/lib:/clusterfs/mp/software/cp2k/elpa-2017.05.002_install/lib/

PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/clusterfs/mp/software/cp2k/plumed2.0_install/lib/pkgconfig

conda activate cms
export LANG=en_GB.UTF-8
export PATH
export LD_LIBRARY_PATH
