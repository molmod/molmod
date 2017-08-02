#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in `find * | egrep '(\.py$)|(\.c$)|(\.pyx$)|(\.pxd$)|(\.h$)|(\.rst$)|(^scripts/)'`; do
  echo $file
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
  sed -i -e $'s/^! --$/!--/' ${file}
  sed -i -e $'s/^# --$/#--/' ${file}
  sed -i -e $'s/^\/\/ --$/\/\/--/' ${file}
  ./updateheaders.py ${file}
done
exit 0
