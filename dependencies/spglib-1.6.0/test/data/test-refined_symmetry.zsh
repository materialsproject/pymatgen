#!/bin/zsh

for j in `/bin/ls -d *(/)`;do
  echo $j
  cd $j
  for i in `ls POSCAR-*`;do
    echo -n $i"  "
    ../../check_consistency.py $i
  done
cd ..
done
