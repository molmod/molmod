#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find doc molmod test examples | egrep "(\.rst$)|(\.py$)|(\.c$)|(\.h$)|(\.i$)|(\.pyf$)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done
