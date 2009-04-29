#! /bin/sh
# This is a very simplistic install scipt. Use with care!

if [ -z $1 ] && [ "$1" = "--system" ]; then
  ./uninstall.sh --system
  (cd ext; python setup.py install)
  python setup.py install
else
  ./uninstall.sh
  (cd ext; python setup.py install --home=$HOME)
  python setup.py install --home=$HOME
  echo "Don't forget to add 'export PYTHONPATH=$HOME' to your .bashrc file."
fi
./cleanfiles.sh
(cd ext; ./cleanfiles.sh)

