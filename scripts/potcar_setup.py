#!/usr/bin/env python

import os
import glob
import shutil
import subprocess

pspdir = ""

count = 0
while not os.path.exists(pspdir):
    if count != 0:
        print("Invalid vasp dir!")
    pspdir = raw_input("Please enter full path where the POT_GGA_PAW_PBE, "
                       "etc. subdirs are present. If you obtained the PSPs "
                       "directly from VASP, this should typically be the "
                       "directory that you untar the files to : ")
    print

targetdir = raw_input("Please enter the fullpath of the where you want to "
                      "create your pymatgen resources directory: ")
print

os.makedirs(targetdir)
print("Generating pymatgen resources directory")

for (parent, subdirs, files) in os.walk(pspdir):
    for subdir in subdirs:
        filenames = glob.glob(os.path.join(parent, subdir, "POTCAR*"))
        if len(filenames) > 0:
            basedir = os.path.join(targetdir, os.path.basename(parent))
            if not os.path.exists(basedir):
                os.makedirs(basedir)
            fname = filenames[0]
            dest = os.path.join(basedir, os.path.basename(fname))
            shutil.copy(fname, dest)
            ext = fname.split(".")[-1]
            if ext.upper() in ["Z", "GZ"]:
                subprocess.Popen(["gunzip", dest]).communicate()
            elif ext.upper() in ["BZ2"]:
                subprocess.Popen(["bunzip2", dest]).communicate()
            if subdir == "Osmium":
                subdir = "Os"
            dest = os.path.join(basedir, "POTCAR.{}".format(subdir))
            shutil.move(os.path.join(basedir, "POTCAR"), dest)
            subprocess.Popen(["gzip", dest]).communicate()

print
print("PSP resources directory generated. You should now add the following to "
      "your environment.")
print("export VASP_PSP_DIR={}".format(os.path.abspath(targetdir)))

