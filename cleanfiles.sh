#!/bin/bash
for i in $(find lib test ext | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done

rm -vr test/output
rm -v test/molmod

rm -v MANIFEST
rm -vr dist
rm -vr build
