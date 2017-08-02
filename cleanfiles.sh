#!/bin/bash
for i in $(find molmod data | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -fv ${i}; done
rm -fvr doc/_build/
rm -fvr output
rm -fv MANIFEST
rm -fvr dist
rm -fvr build
rm -fvr doctrees
rm -fv molmod/ext.c
rm -fv molmod/version.py

# Output of the examples
rm -fv data/examples/001_molecules/ibuprofen.xyz
rm -fv data/examples/001_molecules/ibuprofen_com.xyz
