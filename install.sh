#! /bin/bash
# This is a very simplistic install scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  ./uninstall.sh --system
  (cd ext; python setup.py install --prefix=/usr/local)
  python setup.py install --prefix=/usr/local
  ./cleanfiles.sh
  (cd ext; ./cleanfiles.sh)
else
  ./uninstall.sh
  (cd ext; python setup.py install --home=$HOME)
  python setup.py install --home=$HOME
  ./cleanfiles.sh
  (cd ext; ./cleanfiles.sh)
  echo "Don't forget to add 'export PYTHONPATH=\$HOME/lib/python' to your .bashrc file."
fi

