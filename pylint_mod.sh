#!/bin/bash


MODIFIED_FILES=""
for f in `git diff --diff-filter=ACMRTUXB --name-only HEAD v2020.12.3 | grep -E 'pymatgen.*\.(py)$' | sed '/test_/d' | tr '\n' ' '`
do
    if [ -e $f ]; then
        MODIFIED_FILES="$MODIFIED_FILES $f"
    fi
done
# stop the build if there are Python syntax errors or undefined names
pylint `echo $MODIFIED_FILES`
