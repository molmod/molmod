echo Cleaning python code in \'`pwd`\' and subdirectories
for file in `find * | egrep "(\.py$)|(\.f90$)|(^scripts/tr-)"`; do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
exit 0
