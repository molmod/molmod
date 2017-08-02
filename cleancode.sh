#!/usr/bin/env bash
FILES=$(find . -type f | egrep '(\.py$)|(\.c$)|(\.pyx$)|(\.pxd$)|(\.h$)|(\.rst$)|(\.in$)|(\.yml$)|(\.yaml$)')
for FILE in ${FILES}; do
  echo 'CLEANING' ${FILE}
  sed -i -e $'s/\t/    /' ${FILE}
  sed -i -e $'s/[ \t]\+$//' ${FILE}
  sed -i -e $'s/^! --$/!--/' ${FILE}
  sed -i -e $'s/^# --$/#--/' ${FILE}
  sed -i -e $'s/^\/\/ --$/\/\/--/' ${FILE}
  sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' ${FILE}
done
FILES=$(find . -type f | grep -v molmod/examples | egrep '(\.py$)|(\.c$)|(\.pyx$)|(\.pxd$)|(\.h$)|(\.rst$)|(\.in$)|(\.yml$)|(\.yaml$)')
./updateheaders.py ${FILES}
exit 0
