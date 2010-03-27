#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find doc lib test ext | egrep "(\.rst$)|(\.py$)|(\.c$)|(\.h$)|(\.i$)|(\.pyf$)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
exit 0
