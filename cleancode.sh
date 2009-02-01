#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in `find | egrep "(\.py$)|(\.c$)|(\.h$)|(\.i$)|(\.pyf$)"`; do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
exit 0
