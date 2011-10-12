#! /bin/bash
# This is a very simplistic install scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  ./uninstall.sh --system
  ./cleanfiles.sh
  python setup.py install --prefix=/usr/local
  ./cleanfiles.sh
else
  ./uninstall.sh
  ./cleanfiles.sh
  python setup.py install --home=$HOME
  ./cleanfiles.sh
fi

