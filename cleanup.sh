for i in `find * | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$"` ; do rm -v ${i}; done

rm -vr extensions/build
rm -v extensions/*.so
rm -v extensions/*.c

rm -v MANIFEST
rm -vr dist
rm -vr build
