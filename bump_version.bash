#!/bin/bash 

NEW_VERSION=$1
source VERSION

if [ -z "$NEW_VERSION" ]; then
    echo "a flag for the new version is needed, e.g.: bump_version.bash v1.0.3"
    exit
fi

for file in VERSION README.md src/main.f90; do 
    sed "s/$VERSION/$NEW_VERSION/g" $file > tmp; 
    mv tmp $file
done
