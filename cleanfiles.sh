for i in `find * | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$"` ; do rm -v ${i}; done

rm -vr debian/python-*
rm -vr debian/pycompat
rm -vr debian/compat
rm -vr debian/files
rm -vr debian/stamp-makefile-build
rm -vr python-build-stamp-* 

rm -vr test/output
rm -v test/molmod

rm -v MANIFEST
rm -vr dist
rm -vr build
