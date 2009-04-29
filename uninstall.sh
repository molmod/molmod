#! /bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -vr /usr/share/molmod
  rm -vr /usr/lib/python*/site-packages/molmod
else
  if [ -z $PYTHONPATH ]; then
    echo 'WARNING: $PYTHONPATH is not defined, defaulting to \$HOME/lib/python'
    PYTHONPATH=$HOME/lib/python
  fi
  rm -vr $HOME/share/molmod
  rm -vr $PYTHONPATH/molmod
fi
