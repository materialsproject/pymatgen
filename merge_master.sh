#!/bin/bash

for branch in mapi2 pourbaix
do
    git checkout $branch
    git merge master -m "auto merge"
done
git checkout master
git push
