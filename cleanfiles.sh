#!/bin/bash
for i in $(find molmod test examples | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
rm -vr doc/_build/
rm -vr output
rm -v MANIFEST
rm -vr dist
rm -vr build
rm -vr doctrees
rm -v molmod/extmodule.c

# Output of the examples
rm -v data/examples/001_molecules/ibuprofen.xyz
rm -v data/examples/001_molecules/ibuprofen_com.xyz
