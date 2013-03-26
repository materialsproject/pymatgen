#!/bin/bash

for branch in mapi2 pourbaix
do
    git checkout $branch
    git merge master -m "auto merge"
    git push
done
git checkout master
