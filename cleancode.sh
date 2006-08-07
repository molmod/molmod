echo Cleaning python code in \'`pwd`\' and subdirectories
for file in `find * | grep "\.py$"`; do
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
exit 0
