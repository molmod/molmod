#!/bin/bash
for i in $(find molmod test ext | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done

rm -vr doc/_build/

rm -vr test/output
rm -v test/molmod

rm -v MANIFEST
rm -vr dist
rm -vr build
