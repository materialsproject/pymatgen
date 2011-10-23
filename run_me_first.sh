#!/bin/bash
CONFIG_FILE=pymatgen/pymatgen.cfg

echo "Please enter full path where the POT_GGA_PAW_PBE, etc. subdirs are present. (Leave blank to skip this step):"
read pspdir
if [ "$pspdir" != "" ]; then

    echo "Please enter the fullpath of the where you want to create your pymatgen"
    echo "resources directory."
    read targetdir

    mkdir $targetdir
    echo
    echo "Generating pymatgen resources directory"
    for subdir in `ls -d $pspdir/POT*`
    do
        basedir=`basename $subdir`
        currdir=$targetdir/$basedir
        mkdir $currdir
        
        for pdir in `ls -d $subdir/*`
        do
            pbasedir=`basename $pdir`
            if [ -e "$pdir/POTCAR.Z" ] || [ -e "$pdir/POTCAR.gz" ]; then
                potcarfile=`ls -d $pdir/POTCAR*`
                ext=${potcarfile##*.}
                cp $potcarfile $currdir/POTCAR.$pbasedir.Z
                gunzip $currdir/POTCAR.$pbasedir.Z
            elif [ -e "$pdir/POTCAR.bz2" ]; then
                potcarfile=`ls -d $pdir/POTCAR*`
                ext=${potcarfile##*.}
                cp $potcarfile $currdir/POTCAR.$pbasedir.bz2
                bunzip2 $currdir/POTCAR.$pbasedir.bz2
            elif [ -e "$pdir/POTCAR" ]; then
                potcarfile=`ls -d $pdir/POTCAR*`
                cp $potcarfile $currdir/POTCAR.$pbasedir
            fi
        done
        gzip $currdir/*
    done
    echo "[VASP]" > $CONFIG_FILE
    echo "pspdir=$targetdir" >> $CONFIG_FILE
else
    echo "Skipping PSP setup.  Note that POTCAR creation will be limited without this step."
fi

